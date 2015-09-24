/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils;

import HistogramUtils.BamMetadataSampler;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * @author Derek.Bickhart
 */
public class WindowPlan {
    // diagram: chromosome -> ordered list of start coords
    private final Map<String, Integer[]> starts = new ConcurrentHashMap<>();
    private final Map<String, Integer[]> ends = new ConcurrentHashMap<>();
    
    public void GenerateWindows(BamMetadataSampler bam){
        // For now, I'm going to pretend that there is no difference among sex chromosomes
        // TODO: Implement sex-specific chromosome counts
        double fullXCov = bam.chrXCov.values().stream()
                .mapToDouble(s -> s)
                .average()
                .getAsDouble();
        
        // TODO: Make window size more dynamic dependent on read length and coverage
        int windowSize = 0;
        if(fullXCov < 6.0d){
            windowSize = 500;
        }else{
            windowSize = 100;
        }
        
        for(String chr : bam.chrOrder){
            int chrlen = bam.chrLens.get(chr);
            int numIdx = (int) Math.ceil(chrlen / (double) windowSize);
            
            Integer[] start = new Integer[numIdx];
            Integer[] end = new Integer[numIdx];
            
            int prevEnd = 0;
            for(int x = 0; x < numIdx; x++){
                start[x] = prevEnd;
                int checkEnd = prevEnd + windowSize - 1;
                if(checkEnd > chrlen){
                    end[x] = chrlen;
                }else{
                    end[x] = checkEnd;
                }
                prevEnd += windowSize; 
            }
            
            starts.put(chr, start);
            ends.put(chr, end);
        }
    }
    
    public Integer[] getStarts(String chr){
        return this.starts.get(chr);
    }
    public Integer[] getEnds(String chr){
        return this.ends.get(chr);
    }
}
