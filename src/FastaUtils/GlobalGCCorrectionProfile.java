/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaUtils;

import HistogramUtils.BamMetadataSampler;
import HistogramUtils.ChrHistogram;
import HistogramUtils.ChrHistogramFactory;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import stats.StdevAvg;

/**
 *
 * @author desktop
 */
public class GlobalGCCorrectionProfile {
    private static final Logger log = Logger.getLogger(GlobalGCCorrectionProfile.class.getName());
    private final double[] correction = new double[1000];
    
    public void CalculationGCCorrectionValues(BamMetadataSampler bam, ChrHistogramFactory rd, GCWindowFactory gc, double averageRD){
        // I've gotta store bin values for every GC bin from each chromosome
        Map<Integer, List<Double>> binvalues = new HashMap<>();
        
        // Set all current correction bins to -1
        for(int i = 0; i < 1000; i++){
            this.correction[i] = -1;
        }
        
        // Now to loop through each chromosome
        for(String chr : bam.chrOrder){
            ChrHistogram chisto = rd.getChrHistogram(chr);
            GCHistogram gchisto = gc.getGCHistogram(chr);
            
            // Load data for processing
            chisto.readTemp();
            gchisto.readTemp();
            log.log(Level.INFO, "Reading RD and GC info for " + chr + " to calculate GC profile");
            
            if(chisto.getNumEntries() != gchisto.getNumEntries()){
                log.log(Level.SEVERE, "Error! Number of windows in RD histogram does not equal windows in GC histogram! Skipping: " + chr);
                chisto.clearData();
                gchisto.clearData();
                continue;
            }
            
            for(int n = 0; n < chisto.getNumEntries(); n++){
                int gcidx = gchisto.getRoundedGC(n);
                if(!binvalues.containsKey(gcidx))
                    binvalues.put(gcidx, new ArrayList<>());
                
                binvalues.get(gcidx).add(chisto.getScore(n));
            }
            
            chisto.clearData();
            gchisto.clearData();
        }
        
        // Average everything
        Map<Integer, Double> averages = binvalues.entrySet().stream()
                .collect(Collectors.toConcurrentMap(s -> s.getKey(), (Map.Entry<Integer, List<Double>> s) ->{
                    return StdevAvg.DoubleAvg(s.getValue());
                }));
        
        // Set the known values in the final array
        averages.entrySet().stream()
                .forEach((s) -> {
                    this.correction[s.getKey()] = averageRD / s.getValue();
                });
        
        // TODO: determine if population-based screening needs me to fill in all elements of the array
    }
    
    public ChrHistogramFactory CorrectGC(Path tmpDir, BamMetadataSampler bam, ChrHistogramFactory rd, GCWindowFactory gc){
        ChrHistogramFactory CorrectedRD = new ChrHistogramFactory();
        for(String chr : bam.chrOrder){
            ChrHistogram chisto = rd.getChrHistogram(chr);
            GCHistogram gchisto = gc.getGCHistogram(chr);
            
            // Load data for processing
            chisto.readTemp();
            gchisto.readTemp();
            log.log(Level.INFO, "Reading RD and GC info for " + chr + " to correct GC bias");
            
            if(chisto.getNumEntries() != gchisto.getNumEntries()){
                log.log(Level.SEVERE, "Error! Number of windows in RD histogram does not equal windows in GC histogram! Not performing correction for: " + chr);
                chisto.clearData();
                gchisto.clearData();
                continue;
            }
            
            for(int n = 0; n < chisto.getNumEntries(); n++){
                int gcidx = gchisto.getRoundedGC(n);
                if(this.correction[gcidx] == -1)
                    continue;
                double score = chisto.getScore(n);
                double correctval = this.correction[gcidx];
                double correctedscore = score * correctval;
                if(Double.isNaN(correctedscore)){
                    correctedscore = 0.0d;
                }
                CorrectedRD.addHistogramData(tmpDir, chr, chisto.getStart(n), chisto.getEnd(n), correctedscore);
            }
            
            log.log(Level.INFO, "Completed correction of " + chr + " and now writing corrected RD histogram to tmp file");
            CorrectedRD.getChrHistogram(chr).writeToTemp();
        }
        return CorrectedRD;
    }
}
