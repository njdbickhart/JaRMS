/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Profiler;

import HistogramUtils.ChrHistogramFactory;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author dbickhart
 */
public class SampleList {
    private static final Logger log = Logger.getLogger(SampleList.class.getName());
    private Map<String, ChrHistogramFactory> data;
    private final List<String> sampleNames;
    private List<Path> InputFiles;
    private String outBase;
    
    public SampleList(List<Path> inputFiles, List<String> sampleNames, String outBase){
        this.InputFiles = inputFiles;        
        this.outBase = outBase;
        this.sampleNames = sampleNames;
    }
    
    public void ProcessFiles(){
        
    }
    
    public Set<String> getSamples(){
        return this.data.keySet();
    }
    
    public double[] getRDValues(String Sample, String chr){
        return this.data.get(Sample)
                .getChrHistogram(chr)
                .retrieveRDBins();
    }
}
