/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

// Collect read midpoint position count arrays for each bam

import DataUtils.WindowPlan;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

// Merge all counts into a separate position count array for the bams
// Calculate the GC percentage per chromosome
// Normalize using similar method to lowess smoother in Alkan
// Do mean shift algorithm with three levels and increasing bin sizes -- this makes the partition
// Segment the genome based on the partition signal -- two-pass segmentation in order to join smaller segments
/**
 *
 * @author Derek.Bickhart
 */
public class ChrHistogramFactory {
    private final Map<String, ChrHistogram> histograms = new ConcurrentHashMap<>();
    
    public void processBamNoRG(String BamFile, WindowPlan wins, String tmpdir) throws Exception{
        Path BamPath = Paths.get(BamFile);
        Path tmp = Paths.get(tmpdir);
        
        SamReader reader = SamReaderFactory.make()
                .validationStringency(ValidationStringency.LENIENT)
                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
                .open(BamPath.toFile());
        
        
        AlignmentIteration(reader, wins, tmp);
    }

    private void AlignmentIteration(SamReader reader, WindowPlan wins, Path tmp) throws Exception {
        Integer[] starts = null; Integer[] ends = null;
        String prevChr = "NA"; int curItr = 0; int count = 0; boolean startL = false;
        // Assuming BAM is sorted
        for(SAMRecord s : reader){
            String chr = s.getReferenceName();
            if(!chr.equals(prevChr) && startL){                
                // cleanup memory
                histograms.get(prevChr).addHistogram(chr, starts[curItr], ends[curItr], (double) count);
                starts = wins.getStarts(chr);
                ends = wins.getEnds(chr);
                histograms.get(prevChr).writeToTemp();
                prevChr = chr;
                count = 0;
                histograms.put(chr, new ChrHistogram(chr, tmp));
            }else if(!startL){
                prevChr = chr;
                starts = wins.getStarts(chr);
                ends = wins.getEnds(chr);
                histograms.put(chr, new ChrHistogram(chr, tmp));
                startL = true;
            }
            
            // TODO: Add in better logic for removing supplementary alignments and/or low MQ reads
            int mid = (s.getAlignmentEnd() + s.getAlignmentStart()) / 2;
            if(mid > ends[curItr]){
                histograms.get(chr).addHistogram(chr, starts[curItr], ends[curItr], (double) count);
                count = 0;
                while(mid > ends[curItr]){
                    curItr++;
                    if(curItr > ends.length)
                        throw new Exception("Error! Bam file alignment inconsistent with Chromosome length: " + chr);
                    else if (mid <= ends[curItr] && mid >= starts[curItr])
                        break;
                    // Add empty window -- the bam is sorted and there are no alignments in this window!
                    histograms.get(chr).addHistogram(chr, starts[curItr], ends[curItr], (double) count);
                }
            }
            count++;
        }
        if(curItr > 0){
            // Add final data to histogram if we've already loaded some data previously
            histograms.get(prevChr).addHistogram(prevChr, starts[curItr], ends[curItr], (double) count);
            histograms.get(prevChr).writeToTemp();
        }
    }
    
    public void addHistogramData(Path tmpDir, String chr, int start, int end, double score){
        if(!this.histograms.containsKey(chr)){
            this.histograms.put(chr, new ChrHistogram(chr, tmpDir));
        }
        this.histograms.get(chr).addHistogram(chr, start, end, score);
    }
    
    public ChrHistogram getChrHistogram(String chr){
        return this.histograms.get(chr);
    }
}
