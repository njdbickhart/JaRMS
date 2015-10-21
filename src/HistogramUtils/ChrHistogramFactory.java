/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

// Collect read midpoint position count arrays for each bam

import DataUtils.ThreadingUtils.ThreadHistoFactory;
import DataUtils.ThreadingUtils.ThreadTempRandAccessFile;
import DataUtils.WindowPlan;
import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

// Merge all counts into a separate position count array for the bams
// Calculate the GC percentage per chromosome
// Normalize using similar method to lowess smoother in Alkan
// Do mean shift algorithm with three levels and increasing bin sizes -- this makes the partition
// Segment the genome based on the partition signal -- two-pass segmentation in order to join smaller segments
/**
 *
 * @author Derek.Bickhart
 */
public class ChrHistogramFactory implements ThreadHistoFactory{
    private static final Logger log = Logger.getLogger(ChrHistogramFactory.class.getName());
    private final Map<String, ChrHistogram> histograms = new ConcurrentHashMap<>();
    
    private SamReader reader;
    private WindowPlan wins;
    private ThreadTempRandAccessFile tmp;
    private Set<String> designatedChrs;
    
    public ChrHistogramFactory(ThreadTempRandAccessFile tmp){
        this.tmp = tmp;
    }
    
    public ChrHistogramFactory(SamReader reader, WindowPlan wins, ThreadTempRandAccessFile tmp){
        this.reader = reader;
        this.wins = wins;
        this.tmp = tmp;
    }
    
    public void processBamNoRG(String BamFile, WindowPlan wins, String tmpdir) throws Exception{
        Path BamPath = Paths.get(BamFile);
        Path tmp = Paths.get(tmpdir);
        
        HTSAlignmentIteration(BamPath, wins);
    }
    
    private void HTSAlignmentIteration(Path BamPath, WindowPlan wins) throws IOException{
        
        for(String chr : wins.getChrList()){
            log.log(Level.INFO, "Getting RD alignments for chr: " + chr);
            Integer[] starts = wins.getStarts(chr);
            Integer[] ends = wins.getEnds(chr);
            final double[] score = new double[starts.length];
            for(int i = 0; i < score.length; i++){
                score[i] = 0.0d;
            }
            SamReader reader = SamReaderFactory.make()
                .validationStringency(ValidationStringency.LENIENT)
                .samRecordFactory(DefaultSAMRecordFactory.getInstance())
                .open(BamPath.toFile());
            
            reader.query(chr, starts[0] + 1, ends[ends.length - 1], true)
                    .forEachRemaining(s -> {
                        int mid = (s.getAlignmentEnd() + s.getAlignmentStart()) / 2;
                        int bin = Math.floorDiv(mid, wins.getWindowSize());
                        
                        if(bin < score.length)
                            score[bin] += 1.0d;
                    });
            
            reader.close();
            this.histograms.put(chr, new ChrHistogram(chr, this.tmp));
            for(int i = 0; i < score.length; i++){
                this.histograms.get(chr).addHistogram(chr, starts[i], ends[i], score[i]);
            }
            this.histograms.get(chr).writeToTemp();
        }
    }

    private void AlignmentIteration(SamReader reader, WindowPlan wins, Path tmp) throws Exception {
        ThreadTempRandAccessFile rand = new ThreadTempRandAccessFile(Paths.get(tmp.toString() + ".rdhisto.tmp"));        
        Integer[] starts = null; Integer[] ends = null;
        String prevChr = "NA"; int curItr = 0; int count = 0; boolean startL = false;
        // Assuming BAM is sorted
        for(SAMRecord s : reader){
            String chr = s.getReferenceName();
            if(!chr.equals(prevChr) && startL){                
                // cleanup memory
                if(curItr != -1){
                    histograms.get(prevChr).addHistogram(chr, starts[curItr], ends[curItr], (double) count);
                    curItr++;
                    while(curItr < starts.length){
                        // Account for empty bins at the end of the chromosome
                        histograms.get(prevChr).addHistogram(chr, starts[curItr], ends[curItr], (double) 0.0d);
                        curItr++;
                    }
                }
                starts = wins.getStarts(chr);
                ends = wins.getEnds(chr);
                log.log(Level.FINEST, "Writing chr: " + prevChr + " to temp file...");
                histograms.get(prevChr).writeToTemp();
                prevChr = chr;
                count = 0; curItr = 0;
                histograms.put(chr, new ChrHistogram(chr, rand));
            }else if(!startL){
                prevChr = chr;
                starts = wins.getStarts(chr);
                ends = wins.getEnds(chr);
                histograms.put(chr, new ChrHistogram(chr, rand));
                startL = true;
            }
            
            // TODO: Add in better logic for removing supplementary alignments and/or low MQ reads
            int mid = (s.getAlignmentEnd() + s.getAlignmentStart()) / 2;
            if(curItr == -1){
                continue;
            }
            if(mid > ends[curItr]){
                histograms.get(chr).addHistogram(chr, starts[curItr], ends[curItr], (double) count);
                count = 0;
                while(mid > ends[curItr]){
                    curItr++;
                    if(curItr >= ends.length){
                        log.log(Level.FINEST, "Clipped end of chr: " + chr + "bin: " + curItr);
                        curItr = -1;
                        break;
                    }else if (mid <= ends[curItr] && mid >= starts[curItr])
                        break;
                    // Add empty window -- the bam is sorted and there are no alignments in this window!
                    histograms.get(chr).addHistogram(chr, starts[curItr], ends[curItr], (double) count);
                }
            }
            count++;
        }
        if(curItr > 0){
            // Add final data to histogram if we've already loaded some data previously
            log.log(Level.FINEST, "Writing last chr: " + prevChr + " to temp file...");
            histograms.get(prevChr).addHistogram(prevChr, starts[curItr], ends[curItr], (double) count);
            curItr++;
            while(curItr < starts.length){
                // Account for empty bins at the end of the chromosome
                histograms.get(prevChr).addHistogram(prevChr, starts[curItr], ends[curItr], (double) 0.0d);
                curItr++;
            }
            histograms.get(prevChr).writeToTemp();
        }
    }
    
    public void addHistogramData(ThreadTempRandAccessFile rand, String chr, int start, int end, double score){
        if(!this.histograms.containsKey(chr)){
            this.histograms.put(chr, new ChrHistogram(chr, rand));
        }
        this.histograms.get(chr).addHistogram(chr, start, end, score);
    }
    public boolean hasChrHistogram(String chr){
        return this.histograms.containsKey(chr);
    }
    public ChrHistogram getChrHistogram(String chr){
        return this.histograms.get(chr);
    }
    
    public void checkSumScores(){
        for(String chr : this.histograms.keySet()){
            this.histograms.get(chr).recalculateSumScore();
        }
    }

    @Override
    public void ProcessSpecificWorkload(Set<String> chrs) {
        this.designatedChrs = chrs;
    }

    @Override
    public void Consolidate(ThreadTempRandAccessFile rand) {
        for(String chr : this.histograms.keySet()){
            this.histograms.get(chr).TransferToSharedTemp(rand);
        }
        this.tmp.Close();
    }

    @Override
    public void run() {
        for(String chr : this.designatedChrs){
            Integer[] starts = this.wins.getStarts(chr);
            Integer[] ends = this.wins.getEnds(chr);
            final double[] score = new double[starts.length];
            for(int i = 0; i < score.length; i++){
                score[i] = 0.0d;
            }
            
            this.reader.query(chr, starts[0], ends[ends.length - 1], true)
                    .forEachRemaining(s -> {
                        int mid = (s.getAlignmentEnd() + s.getAlignmentStart()) / 2;
                        int bin = Math.floorDiv(mid, this.wins.getWindowSize());
                        
                        if(bin < score.length)
                            score[bin] += 1.0d;
                    });
            
            this.histograms.put(chr, new ChrHistogram(chr, this.tmp));
            for(int i = 0; i < score.length; i++){
                this.histograms.get(chr).addHistogram(chr, starts[i], ends[i], score[i]);
            }
            this.histograms.get(chr).writeToTemp();
        }
    }

    @Override
    public void ResumeFromTempFile(ThreadTempRandAccessFile rand) {
        for(String chr : this.tmp.getListChrs()){
            this.histograms.put(chr, new ChrHistogram(chr, this.tmp));
            this.histograms.get(chr).setNumEntries((int) (this.tmp.getChrLength(chr) / 16));
            this.histograms.get(chr).recalculateSumScore();
        }
    }
    
    public void PrintWindowsFromTempFile(BufferedWriter output){
        for(String chr : this.tmp.getListChrs()){
            if(this.histograms.containsKey(chr))
                try {
                    this.histograms.get(chr).WriteOutText(output);
                } catch (IOException ex) {
                    log.log(Level.SEVERE, "Error writing text output for chr: " + chr, ex);
                }
        }
    }
}
