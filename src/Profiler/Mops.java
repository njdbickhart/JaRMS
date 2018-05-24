/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Profiler;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Logger;

/**
 *
 * @author dbickhart
 */
public class Mops {
    private static final Logger log = Logger.getLogger(Mops.class.getName());
    
    private final SampleList list;
    // I is the expected, normalized fold change from CN2. This can be modified for binning non-diploid species as well
    private double I[] = {0.025,0.5,1,1.5,2,2.5,3,3.5,4.0};

    private int cycles = 20;
    
    
    public Mops(SampleList list){
        this.list = list;
        
    }
    
    public void EstimateSINI(){
        // Calculate I/NI values for all chromosomes
    }
    
    public void SegmentSINI(){
        // Segment the I/NI calls and form into discrete classes
        // The goal here is to take sINI calls and push them into separate containers -- one for misassemblies and the other for WSSD
        // INI negative windows are then noted and saved for incorporation into downstream pipelines
    }
    
    
}
