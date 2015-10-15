/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils.ThreadingUtils;

import DataUtils.ThreadingUtils.TempHistogram;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Class that gives UCSC binning properties to Histogram classes for fast index retrieval
 * @author Derek.Bickhart
 */
public class HistogramIndexer {
    // Key reference is: chr -> bin integer -> list of idx integers
    private final Map<String, Map<Integer, List<Integer>>> idxMap = new ConcurrentHashMap<>();
    
    private final static int binOffsetsExtended[] = {4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};
    private final static int _binFirstShift = 17;
    private final static int _binNextShift = 3;
    private final static int _binOffsetOldToExtended = 4681;
    
    public void IndexHistogram(TempHistogram t){
        int idx = 0;
        while((idx = t.iterateAdvance()) != -1){
            int bin = this.getBin(t.getStart(), t.getEnd());
            if(!this.idxMap.containsKey(t.chr))
                this.idxMap.put(t.chr, new ConcurrentHashMap<>());
            
            if(!this.idxMap.get(t.chr).containsKey(bin))
                this.idxMap.get(t.chr).put(bin, new ArrayList<>());
            
            this.idxMap.get(t.chr).get(bin).add(idx);
        }
    }
    
    /*
        From James Kent's source code
    */
    private int getBin(int start, int end){
        int startBin = start, endBin = end-1, i;
        startBin >>= _binFirstShift;
        endBin >>= _binFirstShift;
        for (i=0; i< binOffsetsExtended.length; ++i){
            if (startBin == endBin)
                return _binOffsetOldToExtended + binOffsetsExtended[i] + startBin;
            startBin >>= _binNextShift;
            endBin >>= _binNextShift;
        }        
        return 0;
    }
}
