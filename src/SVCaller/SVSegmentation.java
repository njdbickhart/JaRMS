/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SVCaller;

import DataUtils.BedStats;
import DataUtils.BinCoords;
import DataUtils.CallEnum;
import DataUtils.WindowPlan;
import HistogramUtils.ChrHistogramFactory;
import HistogramUtils.LevelHistogram;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import stats.EvalueTools;
import stats.StdevAvg;
import stats.TDistributionFunction;

/**
 *
 * @author desktop
 */
public class SVSegmentation {
    private static final Logger log = Logger.getLogger(SVSegmentation.class.getName());
    private List<BedStats> StoredCalls;
    
    private final Path outfile;
    //private final List<BedStats> CNVCalls = new ArrayList<>();
    
    public SVSegmentation(String outDir){
        this.outfile = Paths.get(outDir + "cnvs.bed");
    }
    
    public void RunSegmentation(final ChrHistogramFactory corrGC, WindowPlan wins, MeanShiftMethod meanShift, int threads){
        // Calculate global mean and SD
        double globalSum = wins.getChrList().stream()
                .map(s -> corrGC.getChrHistogram(s).getSum())
                .reduce(0.0d, (a, b) -> a + b);
        int globalNumBins = wins.getChrList().stream()
                .map(s -> corrGC.getChrHistogram(s).getNumEntries())
                .reduce(0, (a, b) -> a + b);
        
        final double globalMean = globalSum / (double) globalNumBins;
        
        double globalSS = wins.getChrList().stream()
                .map(s -> corrGC.getChrHistogram(s).getSumSquares(globalMean))
                .reduce(0.0d, (a, b) -> a + b);
        globalSS /= globalNumBins - 1;
        
        final double globalStdev = Math.sqrt(globalSS);
        
        // Process chromosomes
        List<Future<List<BedStats>>> SegFutures = new ArrayList<>();
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        TDistributionFunction ttest = new TDistributionFunction(wins.getWindowSize());
        
        for(String chr : wins.getChrList()){
            SegFutures.add(executor.submit(new Segmentation(corrGC.getChrHistogram(chr).retrieveRDBins(), 
                    meanShift.getChrLevels(chr), 
                    ttest, 
                    wins, 
                    chr, 
                    globalMean, 
                    globalStdev)));
        }
        
        executor.shutdown();
        
        while(!executor.isTerminated()){}
        
        List<BedStats> FinalCalls = new ArrayList<>();
        for(Future<List<BedStats>> threadSegs : SegFutures){
            try {
                FinalCalls.addAll(threadSegs.get());
            } catch (InterruptedException | ExecutionException ex) {
                log.log(Level.SEVERE, "Error retrieving Threaded results!", ex);
            }
        }
        
        Collections.sort(FinalCalls);
        this.StoredCalls = FinalCalls;
    }
    
    public void printOutAllCalls(){
        try(BufferedWriter out = Files.newBufferedWriter(outfile, Charset.defaultCharset())){
            for(BedStats b : this.StoredCalls){
                out.write(b.getOutputString());
            }
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error printing out CNV calls!", ex);
        }
    }
    
    
    private class Segmentation implements Callable<List<BedStats>>{
        private final double[] corrRD;
        private final LevelHistogram level;
        private final WindowPlan wins;
        private final double globalMean;
        private final double globalSD;
        private final String chr;
        private final List<CallEnum> calls;
        private final TDistributionFunction ttest;
        
        /*
            Constants
        */
        private final int GenomeSize;
        private final double CUTOFF_REGION = 0.05d;
        private final double CUTOFF_TWO_REGIONS = 0.01;
        private final double GENOME_SIZE_CNV;
        private final EnumSet<CallEnum> cnvEnums = EnumSet.of(CallEnum.DELETION, CallEnum.DUPLICATION, CallEnum.POTENTIALDEL);
        
        private double chrMean;
        private double chrSD;
        
        public Segmentation(double[] corrRD, LevelHistogram level, TDistributionFunction ttest, WindowPlan wins, String chr, double globalMean, double globalSD){
            this.corrRD = corrRD;
            this.level = level;
            this.wins = wins;
            this.chr = chr;
            this.globalMean = globalMean;
            this.globalSD = globalSD;
            this.calls = new ArrayList<>(corrRD.length);
            this.ttest = ttest;
            this.GenomeSize = wins.getGenomeSize();
            this.GENOME_SIZE_CNV = this.GenomeSize * 0.1d;
        }
        
        @Override
        public List<BedStats> call() throws Exception {
            List<BedStats> FinalCalls = new ArrayList<>();
            // calculate mean and SD
            this.chrMean = StdevAvg.DoubleAvg(corrRD);
            this.chrSD = StdevAvg.stdevDBL(chrMean, corrRD);
            
            // merge level signal
            double cut = this.chrMean / 4;
            this.level.performMerger(cut);
            
            // prepare call list
            for(int i = 0; i < this.corrRD.length; i++)
                this.calls.add(CallEnum.NORMAL);
            
            // Go through routines to fill out calls List
            this.IdentifyInitialRegions(cut);
            this.MergeSmallRegions();
            this.IdentifyPotentialDels(this.chrMean - cut);
            
            // Prepare final calls
            for (int b = 0;b < this.corrRD.length; b++) {
                CallEnum c = this.calls.get(b);
                
                if (c.equals(CallEnum.NORMAL)) continue;
                int bs = b;
                double cnv = 0;
                while (b < this.corrRD.length && this.calls.get(b).equals(c)) 
                    cnv += this.corrRD[b++];
                int be = --b;

                if (be <= bs) continue;

                cnv /= (be - bs + 1) * this.chrMean;
                String type;
                switch(c){
                    case DELETION:
                    case POTENTIALDEL:
                        type = "deletion";
                        break;
                    case DUPLICATION:
                        type = "duplication";
                        break;
                    default:
                        type = "???";
                }
                int start  = bs * wins.getWindowSize() + 1;
                int end    = (be + 1) * wins.getWindowSize();
                double e  = EvalueTools.GetEValue(this.chrMean,this.chrSD,this.GenomeSize, this.ttest, this.corrRD,bs,be);
                double e2 = EvalueTools.GetGaussianEValue(this.chrMean,this.chrSD,this.corrRD,bs,be, this.GenomeSize);
                double e3 = 1,e4 = 1;
                int add = Math.round(1000.0f / wins.getWindowSize() + 0.5f);
                if (bs + add < be - add) {
                    //e3 = getEValue(mean,sigma,rd,bs + add,be - add);
                    //e4 = gaussianEValue(mean,sigma,rd,bs + add,be - add);
                }
                double n_reads_all = 0,n_reads_unique = 0;
                for (int i = bs;i <= be;i++) {
                    //n_reads_all    += h_all->GetBinContent(i);
                    //n_reads_unique += h_unique->GetBinContent(i);
                }
                double q0 = -1;
                if (n_reads_all > 0) q0 = (n_reads_all - n_reads_unique)/n_reads_all;
                FinalCalls.add(new BedStats(this.chr, start, end, type, cnv, e, e2));
            }
            
            return FinalCalls;
        }
        
        private void IdentifyInitialRegions(double cut){
            double min =  this.chrMean - cut;
            double max = this.chrMean + cut;
            for (int b = 0;b < this.corrRD.length;b++) {
                BinCoords bin = new BinCoords();
                int b0 = b;
                bin.start = b;
                while (b < this.corrRD.length && level.getScore(b) < min) 
                    b++;
                bin.end = b - 1;
                if (bin.isNormalInterval() && EvalueTools.AdjustToEvalue(chrMean,chrSD, GenomeSize, ttest, this.corrRD,bin,CUTOFF_REGION))
                    for (int i = bin.start;i <= bin.end;i++) 
                        this.calls.set(i, CallEnum.DELETION);
                bin.start = b;
                while (b < this.corrRD.length && level.getScore(b) > max) 
                    b++;
                bin.end = b - 1;
                if (bin.isNormalInterval() && EvalueTools.AdjustToEvalue(chrMean,chrSD,GenomeSize, ttest, this.corrRD, bin,CUTOFF_REGION))
                    for (int i = bin.start;i <= bin.end;i++) 
                        this.calls.set(i, CallEnum.DUPLICATION);
                if (b > b0) 
                    b--;
            }
        }
        
        private void MergeSmallRegions(){
            int n_add = 1;
            while (n_add > 0) {
                for (int b = 0;b < this.corrRD.length; b++) {
                    if (this.calls.get(b).equals(CallEnum.NORMAL)) 
                        continue;

                    int s = b;
                    while (b < this.corrRD.length && this.calls.get(b).equals(CallEnum.NORMAL)) 
                        b++;
                    int e = b - 1;

                    if (e < s || s == 0 || e >= this.corrRD.length) 
                        continue;
                    if(this.calls.get(s - 1).equals(this.calls.get(e + 1)))
                        continue;
                    if (s == e) { 
                        this.calls.set(s, this.calls.get(s - 1));
                        continue; 
                    }

                    int le = s - 1,ls = le;
                    while( ls >= 0 && this.calls.get(ls).equals(this.calls.get(le)))
                        ls--;
                    ls++;
                    int rs = e + 1,re = rs;
                    while(re < this.corrRD.length && this.calls.get(re).equals(this.calls.get(rs)))
                        re++;
                    re--;

                    double average = StdevAvg.getRangeAverage(corrRD, s, e), 
                            variance = StdevAvg.getRangeVariance(corrRD, average, s, e);
                    double raverage = StdevAvg.getRangeAverage(corrRD, rs, re), 
                            rvariance = StdevAvg.getRangeVariance(corrRD, raverage, rs, re);
                    double laverage = StdevAvg.getRangeAverage(corrRD, ls, le),
                            lvariance = StdevAvg.getRangeVariance(corrRD, laverage, ls, le);
                    int n = e - s,rn = re - rs, ln = le - ls;
                    if (n > rn || n > ln) 
                        continue;

                    if (this.ttest.TestTwoRegions(laverage,lvariance,ln,average,variance,n,
                                       GENOME_SIZE_CNV) < CUTOFF_TWO_REGIONS &&
                        this.ttest.TestTwoRegions(raverage,rvariance,rn,average,variance,n,
                                       GENOME_SIZE_CNV) < CUTOFF_TWO_REGIONS)
                        continue;
                    
                    for (int i = s; i <= e; i++)
                        this.calls.set(i, CallEnum.POTENTIAL);
                }

                n_add = 0;
                for (int i = 0;i <= this.corrRD.length; i++){ 
                    if(this.calls.get(i).equals(CallEnum.POTENTIAL)){
                        this.calls.set(i, this.calls.get(i -1));
                        n_add++;
                    }
                }
            }
        }
        
        private void IdentifyPotentialDels(double min){
            for (int b = 0;b < this.corrRD.length; b++) {
                if(this.calls.get(b).equals(CallEnum.NORMAL))
                    continue;
                int bs = b;
                while (b < this.corrRD.length && level.getScore(b) < min) b++;
                int be = b - 1;
                if (be > bs) {
                    if (EvalueTools.GetGaussianEValue(this.chrMean,this.chrSD,this.corrRD,bs,be, this.GenomeSize) < CUTOFF_REGION)
                        for(int i = bs; i <= be; i++)
                            this.calls.set(i, CallEnum.POTENTIALDEL);
                    b--;
                }
            }
        }
    }
}
