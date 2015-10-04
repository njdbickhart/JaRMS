/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexMetaData;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 *
 * @author Derek.Bickhart
 */
public class BamMetadataSampler {
    public Map<String, Integer> chrLens;
    public Map<String, Double> chrXCov;
    public Set<String> chrOrder;
    private final Path bamFile;
    
    public BamMetadataSampler(String bamFile){
        this.bamFile = Paths.get(bamFile);
        //if(!this.bamFile.toFile().canRead())
            //throw new Exception("Error! Cannot find BAM file: " + bamFile);
    }
    
    public void getMetaData(){
        SamReader reader = SamReaderFactory.make()
                .validationStringency(ValidationStringency.LENIENT)
                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
                .open(this.bamFile.toFile());
        
        final SAMFileHeader head = reader.getFileHeader();
        
        this.chrLens = head.getSequenceDictionary().getSequences().stream()
                .collect(Collectors.toMap(s -> s.getSequenceName(), s -> s.getSequenceLength()));
        
        this.chrOrder = head.getSequenceDictionary().getSequences().stream()
                .sequential()
                .map(s -> s.getSequenceName())
                .collect(Collectors.toSet());
        
        final BAMIndex index = reader.indexing().getIndex();
        
        this.chrXCov = chrLens.entrySet().stream()
                .collect(Collectors.toMap(s -> s.getKey(), 
                        (Map.Entry<String, Integer> s) -> { 
                            BAMIndexMetaData meta = index.getMetaData(head.getSequenceIndex(s.getKey()));
                            return (meta.getAlignedRecordCount() + meta.getUnalignedRecordCount()) / (double) s.getValue();
                        }));
    }
}
