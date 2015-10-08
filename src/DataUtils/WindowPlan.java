/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils;

import HistogramUtils.BamMetadataSampler;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * @author Derek.Bickhart
 */
public class WindowPlan {
    // diagram: chromosome -> ordered list of start coords
    private final Map<String, Integer[]> starts = new ConcurrentHashMap<>();
    private final Map<String, Integer[]> ends = new ConcurrentHashMap<>();
    private final Map<String, Integer> numBins = new ConcurrentHashMap<>();
    private int windowSize;
    private long GenomeSize;
    
    public void GenerateWindows(BamMetadataSampler bam){
        this.GenomeSize = bam.chrLens.values().stream()
                .reduce(0, (a, b) -> a + b);
        
        // For now, I'm going to pretend that there is no difference among sex chromosomes
        // TODO: Implement sex-specific chromosome counts
        double fullXCov = bam.chrXCov.values().stream()
                .mapToDouble(s -> s)
                .average()
                .getAsDouble();
        
        // TODO: Make window size more dynamic dependent on read length and coverage
        windowSize = 0;
        if(fullXCov < 6.0d){
            windowSize = 500;
        }else{
            windowSize = 100;
        }
        
        bam.chrOrder.stream().forEach((chr) -> {
            int binNums = 0;
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
                binNums++;
            }
            
            starts.put(chr, start);
            ends.put(chr, end);
            binNums++;
            
            this.numBins.put(chr, binNums);
        });
    }
    
    public int getBinStart(String chr, int binNum) throws Exception{
        if(this.numBins.get(chr) > binNum){
            return this.starts.get(chr)[binNum];
        }else{
            throw new Exception("Error! Attempted to retrieve start coordinate outside of bounds! binnum: " + binNum + " Max bin: " + this.numBins.get(chr));
        }
    }
    public int getBinEnd(String chr, int binNum) throws Exception{
        if(this.numBins.get(chr) > binNum){
            return this.ends.get(chr)[binNum];
        }else{
            throw new Exception("Error! Attempted to retrieve end coordinate outside of bounds! binnum: " + binNum + " Max bin: " + this.numBins.get(chr));
        }
    }
    
    public Integer[] getStarts(String chr){
        return this.starts.get(chr);
    }
    public Integer[] getEnds(String chr){
        return this.ends.get(chr);
    }
    public int getWindowSize(){
        return this.windowSize;
    }
    public long getGenomeSize(){
        return this.GenomeSize;
    }
    public Set<String> getChrList(){
        return this.starts.keySet();
    }
    public boolean containsChr(String chr){
        return this.starts.containsKey(chr);
    }
}
