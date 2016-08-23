/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

import DataUtils.WindowPlan;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

/**
 *
 * @author Derek.Bickhart
 */
public class BamMetadataSampler {
    private static final Logger log = Logger.getLogger(BamMetadataSampler.class.getName());
    public Map<String, Long> chrLens = null;
    public final Map<String, Double> chrXCov;
    public Set<String> chrOrder = null;
    private final List<Path> bamFiles;
    
    public BamMetadataSampler(List<String> bamFiles){
        this.bamFiles = bamFiles.stream()
                .map((s) -> Paths.get(s))
                .collect(Collectors.toList());
        this.chrXCov = new HashMap<>();
        //if(!this.bamFile.toFile().canRead())
            //throw new Exception("Error! Cannot find BAM file: " + bamFile);
    }
    
    public void getMetaData(){
        double totalXCov = 0.0d;
        long totalGenomeSize = 0l;
        
        for(Path bam : this.bamFiles){
            SamReader reader = SamReaderFactory.make()
                    .validationStringency(ValidationStringency.LENIENT)
                    .samRecordFactory(DefaultSAMRecordFactory.getInstance())
                    .open(bam.toFile());

            final SAMFileHeader head = reader.getFileHeader();

            if(this.chrLens == null)
                this.chrLens = head.getSequenceDictionary().getSequences().stream()
                        .collect(Collectors.toMap(s -> s.getSequenceName(), s -> (long)s.getSequenceLength()));
            
            if(this.chrOrder == null)
                this.chrOrder = head.getSequenceDictionary().getSequences().stream()
                        .sequential()
                        .map(s -> s.getSequenceName())
                        .collect(Collectors.toSet());
            else{
                Set<String> sequences = head.getSequenceDictionary().getSequences().stream()
                        .sequential()
                        .map(s -> s.getSequenceName())
                        .collect(Collectors.toSet());
                
                for(String chr : this.chrOrder){
                    if(!sequences.contains(chr)){
                        log.log(Level.SEVERE, "Error processing Bam file: " + bam + "! Chr: " + chr + " not in previous file!");
                        System.exit(-1);
                    }
                }
            }
                

            final BAMIndex index = reader.indexing().getIndex();

            chrLens.entrySet().stream()
                    .forEach((s) -> {
                        if(!this.chrXCov.containsKey(s.getKey()))
                            this.chrXCov.put(s.getKey(), 0.0d);
                        BAMIndexMetaData meta = index.getMetaData(head.getSequenceIndex(s.getKey()));
                        if(meta.getAlignedRecordCount() == 0)
                            log.log(Level.FINE, "Found no aligned reads in bam: " + bam + " for chr: " + s.getKey());
                        
                        // Note: this is a rolling average
                        // Also, I'm assuming an average read length of 100
                        this.chrXCov.put(s.getKey(), 
                                (this.chrXCov.get(s.getKey()) + ((meta.getAlignedRecordCount() + meta.getUnalignedRecordCount()) * 100 / (double) s.getValue()))
                                    / 2.0d);
                    });

            totalXCov += this.chrXCov.entrySet().stream()
                    .map(s -> s.getValue())
                    .reduce(0.0d, (a, b) -> (a + b));
            
            if(totalGenomeSize == 0l)
                totalGenomeSize = this.chrLens.entrySet().stream()
                        .map(s -> s.getValue())
                        .reduce(0l, (a, b) -> (a + b));
            
        }
        
        // I realize that I've already factored out the genome size in the chr length estimation
        //totalXCov /= totalGenomeSize;
        totalXCov /= chrOrder.size();
        log.log(Level.INFO, "BAM metadata stats: genome size : " + totalGenomeSize + " #chrs: " + this.chrOrder.size() + " avg X coverage: " + totalXCov);
    }
    
    public void UpdateChrOrder(WindowPlan wins){
        int start = this.chrOrder.size();
        this.chrOrder.removeAll(wins.getExclude());
        log.log(Level.INFO, "Updated mapped chrs. Starting count: " + start + " ending count: " + this.chrOrder.size());
    }
}
