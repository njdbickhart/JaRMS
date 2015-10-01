/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package SVCaller;

import DataUtils.BinCoords;
import DataUtils.CallEnum;
import DataUtils.WindowPlan;
import HistogramUtils.LevelHistogram;
import file.BedSimple;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;
import java.util.concurrent.Callable;
import stats.StdevAvg;
import stats.TDistributionFunction;

/**
 *
 * @author desktop
 */
public class SVSegmentation {
    
    private class Segmentation implements Callable<List<BedSimple>>{
        private final double[] corrRD;
        private final LevelHistogram level;
        private final WindowPlan wins;
        private final double globalMean;
        private final double globalSD;
        private final String chr;
        private final List<EnumSet<CallEnum>> calls;
        private final TDistributionFunction ttest;
        private final int GenomeSize;
        private final double CUTOFF_REGION = 0.05d;
        
        private double chrMean;
        private double chrSD;
        
        public Segmentation(double[] corrRD, LevelHistogram level, TDistributionFunction ttest, int GenomeSize, WindowPlan wins, String chr, double globalMean, double globalSD){
            this.corrRD = corrRD;
            this.level = level;
            this.wins = wins;
            this.chr = chr;
            this.globalMean = globalMean;
            this.globalSD = globalSD;
            this.calls = new ArrayList<>(corrRD.length);
            this.ttest = ttest;
            this.GenomeSize = GenomeSize;
        }
        
        @Override
        public List<BedSimple> call() throws Exception {
            List<BedSimple> calls = new ArrayList<>();
            // calculate mean and SD
            this.chrMean = StdevAvg.DoubleAvg(corrRD);
            this.chrSD = StdevAvg.stdevDBL(chrMean, corrRD);
            
            // merge level signal
            double cut = this.chrMean / 4;
            this.level.performMerger(cut);
            
            
            
            return calls;
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
                if (be > bs && adjustToEValue(chrMean,chrSD, GenomeSize, ttest, this.corrRD,n_bins,bs,be,CUTOFF_REGION))
                    for (int i = bs;i <= be;i++) flags[i] = 'D';
                bs = b;
                while (b < n_bins && level[b] > max) b++;
                be = b - 1;
                if (be > bs && adjustToEValue(chrMean,chrSD,rd,n_bins,bs,be,CUTOFF_REGION))
                    for (int i = bs;i <= be;i++) flags[i] = 'A';
                if (b > b0) b--;
            }
        }
        
    }
}
