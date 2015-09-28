/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SVCaller;

import DataUtils.WindowPlan;
import HistogramUtils.ChrHistogram;
import HistogramUtils.ChrHistogramFactory;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;
import stats.StdevAvg;
import stats.TDistributionFunction;

/**
 *
 * @author desktop
 */
public class MeanShiftMethod {
    private final Map<String, ChrHistogram> shiftedChrHistos = new ConcurrentHashMap<>();
    
    public void Partition(ChrHistogramFactory chisto, WindowPlan wins, Path tmpDir, int range, int threads){
        TDistributionFunction ttest = new TDistributionFunction(wins.getWindowSize());
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        List<Future<ChrHistogram>> workers = new ArrayList<>();
        for(String chr : wins.getChrList()){
            workers.add(executor.submit(new MeanShifter(wins, ttest, chisto.getChrHistogram(chr).retrieveRDBins(), tmpDir, chr, range)));
        }
        
        executor.shutdown();
        while(!executor.isTerminated()){
            
        }
        
        for(Future<ChrHistogram> chrHisto : workers){
            try {
                ChrHistogram c = chrHisto.get();
                this.shiftedChrHistos.put(c.getChr(), c);
            } catch (InterruptedException | ExecutionException ex) {
                Logger.getLogger(MeanShiftMethod.class.getName()).log(Level.SEVERE, "Error retrieving ChrHistogram from threaded worker!", ex);
            }
        }
    }
    
    private class MeanShifter implements Callable<ChrHistogram>{
        private final Logger log = Logger.getLogger(MeanShifter.class.getName());
        
        private final WindowPlan wins;
        private final double[] rdBins;
        private final String chr;
        private final Path tmpDir;
        private final TDistributionFunction ttest;
        // Range is the bin_range which is typically 128
        private final int range;
        private final int invnum = 10000;
        private final double CUTOFF_REGION = 0.05;
        private final double CUTOFF_TWO_REGIONS = 0.01;
        private final double[] inversions = new double[invnum];
        
        public MeanShifter(WindowPlan wins, TDistributionFunction ttest, double[] rdBins, Path tmpDir, String chr, int range){
            this.wins = wins;
            this.rdBins = rdBins;
            this.tmpDir = tmpDir;
            this.chr = chr;
            this.range = range;
            this.ttest = ttest;
            for(int i = 0; i < invnum; i++)
                this.inversions[i] = 0.0d;
        }
        
        @Override
        public ChrHistogram call() throws Exception {
            ChrHistogram shifted = new ChrHistogram(this.chr, this.tmpDir);
            boolean[] mask = new boolean[this.rdBins.length];
            for(int i = 0; i < rdBins.length; i++)
                mask[i] = false;
            
            double[] level = Arrays.copyOf(rdBins, rdBins.length);
            double mean = StdevAvg.DoubleAvg(rdBins);
            double sigma = StdevAvg.stdevDBL(mean, rdBins);
            
            for (int bin_band = 2;bin_band <= range;bin_band++) {
      
                log.log(Level.INFO, "Setting bin band for chr: " + chr + " to: " + bin_band);

                for (int b = 0;b < this.rdBins.length; b++) 
                    if (!mask[b]) level[b] = rdBins[b];

                CalcLevels(level, mask, bin_band, mean, sigma);

                CalcLevels(level, mask, bin_band, mean, sigma);

                CalcLevels(level, mask, bin_band, mean, sigma);

                
                UpdateMask(level, mask, mean, sigma);
                

                if (bin_band >=   8) bin_band +=  1;
                if (bin_band >=  16) bin_band +=  2;
                if (bin_band >=  32) bin_band +=  4;
                if (bin_band >=  64) bin_band +=  8;
                if (bin_band >= 128) bin_band += 16;
                if (bin_band >= 256) bin_band += 32;
                if (bin_band >= 512) bin_band += 64;
            }
            // Create chromosome histogram and return
            for(int i = 0; i < level.length; i++){
                shifted.addHistogram(chr, wins.getBinStart(chr, i), wins.getBinEnd(chr, i), level[i]);
            }
            shifted.writeToTemp();
            return shifted;
        }
        
        private void UpdateMask(double[] level, boolean[] mask, double mean, double sigma){
            for (int b = 0;b < this.rdBins.length; b++) mask[b] = false;

            int ln = 0,n = 0,rn = 0;
            double average  = 0, variance  = 0;
            double laverage = 0, lvariance = 0;
            double raverage = 0, rvariance = 0;
            double inv_mean = 1 / mean;
            BinRange curRange = new BinRange();
            while (getRegionRight(level, curRange.stop, curRange)) {

                ln         = n;
                laverage   = average;
                lvariance  = variance;

                n        = rn;
                average  = raverage;
                variance = rvariance;

                BinRange rightRange = new BinRange();
                //int rstart,rstop;
                if (!getRegionRight(level, curRange.stop + 1, rightRange)) break;
                raverage = StdevAvg.getRangeAverage(this.rdBins, rightRange.start, rightRange.stop);
                rvariance = StdevAvg.getRangeVariance(this.rdBins, raverage, rightRange.start, rightRange.stop);

                BinRange leftRange = new BinRange();
                if (!getRegionLeft(level, curRange.start - 1, leftRange)) {
                    // Have no left region -- need to calculate variances
                    average = StdevAvg.getRangeAverage(this.rdBins, curRange.start, curRange.stop);
                    variance = StdevAvg.getRangeVariance(this.rdBins, average, curRange.start, curRange.stop);
                    continue;
                }

                if (n <= 1) 
                    continue;

                if (ln <= 15 || n <= 15 || rn <= 15) {
                    // Check sigma condition
                    double ns = 1.8;
                    //double ns = 2.0;
                    double nsigma = ns * Math.sqrt(level[leftRange.stop] * inv_mean) * sigma;
                    if (Math.abs(level[leftRange.stop] - level[curRange.start]) < nsigma) 
                        continue;
                    nsigma = ns * Math.sqrt(level[rightRange.start] * inv_mean) * sigma;
                    if (Math.abs(level[rightRange.start] - level[curRange.stop]) < nsigma) 
                        continue;
                } else {
                    // Checking compartibility of regions
                    if (ttest.TestTwoRegions(laverage,lvariance,ln,average,variance,n,
                                       wins.getGenomeSize()) > CUTOFF_TWO_REGIONS ||
                        ttest.TestTwoRegions(raverage,rvariance,rn,average,variance,n,
                                       wins.getGenomeSize()) > CUTOFF_TWO_REGIONS)
                      continue;
                }

                // Check that region is abnormal (deviates largerly from the mean)
                if (ttest.TestOneRegion(mean,average,variance,n) > CUTOFF_REGION) 
                    continue;

                // Mask
                for (int b = curRange.start; b <= curRange.stop; b++) 
                    mask[b] = true;
            }
        }
        
        private boolean getRegionRight(double[] level, int bin, BinRange range){
            if (bin < 0 || bin >= this.rdBins.length) 
                return false;
            int start = bin, stop = start;
            while (stop < this.rdBins.length && sameLevel(level[start],level[stop])) 
                stop++;
            stop--;
            range.start = start;
            range.stop = stop;
            return true;
        }
        
        private boolean getRegionLeft(double[] level, int bin, BinRange range){
            if (bin < 0 || bin >= this.rdBins.length) 
                return false;
            int start = bin, stop = bin;
            while (start >= 0 && sameLevel(level[start],level[stop])) 
                start--;
            start++;
            range.start = start;
            range.stop = stop;
            return true;
        }
        
        private boolean sameLevel(double t, double l){
            return Math.abs(t - l) < 0.01d;
        }
        
        private class BinRange{
            public int start = 0;
            public int stop = 0;
        }
        
        private void CalcLevels(double[] level, boolean[] mask, int bin_band, double mean, double sigma){
            double[] grad_b = new double[this.rdBins.length];
            for (int b = 0;b < this.rdBins.length;b++) 
                grad_b[b] = 0;

            double inv2_bin_band = 1.0d/(bin_band * bin_band);
            double mean_4 = mean/4, sigma_2 = 4/(sigma * sigma), ms2 = mean/(sigma * sigma);
            int win = 3*bin_band;
            double exps[] = new double[win + 1];
            for (int i = 0;i <= win;i++)
                exps[i] = i * Math.exp(-0.5 * i * i * inv2_bin_band);
            for (int b = 0;b < this.rdBins.length;b++) {
                if (mask[b]) 
                    continue;
                double inv_b; 
                int d = 0;
                if (level[b] < mean_4) 
                    inv_b = sigma_2;
                else 
                    inv_b = ms2/level[b];
                for (int i = b + 1;i < this.rdBins.length;i++) {
                    if (mask[i]) continue;
                    d++;
                    double inv_i;
                    if (level[i] < mean_4) inv_i = sigma_2;
                    else                   inv_i = ms2/level[i];
                    double r = level[i] - level[b];
                    double val = -0.5*r*r;
                    grad_b[b] += exps[d] * Math.exp(val*inv_b);
                    grad_b[i] -= exps[d] * Math.exp(val*inv_i);
                    if (d == win) break;
                }
            }

            // Calculating levels
            for (int b = 0;b < this.rdBins.length;b++) {
                if (mask[b]) continue;
                int b_start = b;

                // Finding region by bins
              
                while (b < this.rdBins.length && grad_b[b] >= 0 && !mask[b]) b++;
                while (b < this.rdBins.length && grad_b[b] <  0 && !mask[b]) b++;
              
                int b_stop = --b;
                if (b_start > b_stop) {
                    log.log(Level.WARNING, "Abnormal bin range! " + b_start + " is greater than: " + b_stop);
                    b = b_start;
                    continue;
                }

                // Calculating level
                double nl = 0;
                int n = 0;
                for (int i = b_start;i <= b_stop;i++) {
                    if (mask[i]) continue;
                    nl += level[i];
                    n++;
                }
                if (n <= 0) {
                    log.log(Level.WARNING, "Error! Found width of 0!");
                    continue;
                }
                nl *= getInverse(n);
                for (int i = b_start;i <= b_stop;i++) 
                    if (!mask[i]) level[i] = nl;
            }
        }
        
        private double getInverse(int n){
            if (n <= 0) return 0;
            if (n >= this.invnum) return 1./n;
            if (this.inversions[n] == 0) this.inversions[n] = 1./n;
            return this.inversions[n];
        }
    }
}