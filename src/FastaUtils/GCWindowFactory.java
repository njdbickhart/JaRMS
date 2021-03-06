/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaUtils;

import DataUtils.ThreadingUtils.ThreadTempRandAccessFile;
import DataUtils.WindowPlan;
import HistogramUtils.BamMetadataSampler;
import TempFiles.binaryUtils.IntUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * This class estimates GC content in fasta files
 * @author Derek.Bickhart
 */
public class GCWindowFactory {
    private static final Logger log = Logger.getLogger(GCWindowFactory.class.getName());
    
    protected final Map<String, GCHistogram> histograms = new HashMap<>();
    protected final Path fastaPath;
    protected final Path tmpPath;
    protected String expectedFastaProfile;

    public GCWindowFactory(String fastaFile, String tmpdir){
        this.fastaPath = Paths.get(fastaFile);
        this.tmpPath = Paths.get(tmpdir);
    }
    
    public void generateGCProfile(BamMetadataSampler bamMeta, WindowPlan wins){
        // Check if we can load existing profile file
        if(this.checkIfProfileFileExists(wins.getWindowSize())){
            this.generateGCHistoFromProfile(bamMeta);
        }else{
            this.processFastaFile(wins);
            this.writeGCProfileOut(bamMeta);
        }        
    }
    
    private void processFastaFile(WindowPlan wins){
        ThreadTempRandAccessFile rand = new ThreadTempRandAccessFile(Paths.get(this.tmpPath.toString() + ".gcprofile.tmp"));
        try(BufferedReader fasta = Files.newBufferedReader(fastaPath, Charset.defaultCharset())){
            log.log(Level.FINEST, "Beginning fasta file processing");
            String line, prevChr = "NA";
            Integer[] starts = null, ends = null; 
            boolean skip = false;
            int totalCount = 0, bpCount = 0, binIdx = 0, gcCount = 0;
            while((line = fasta.readLine()) != null){
                if(line.startsWith(">")){
                    skip = false;
                    line = line.trim();
                    line = line.replaceAll(">", "");
                    if(!wins.containsChr(line)){
                        log.log(Level.WARNING, "Warning! Bam file either doesn't contain alignments from chr " + line + " or it is too small for analysis!");
                        if(!prevChr.equals("NA")){
                            if(binIdx < starts.length){
                                // Avoids issue where the binIdx is too high because of exact multiple of window size for the fasta entry
                                this.histograms.get(prevChr).addHistogram(prevChr, starts[binIdx], ends[binIdx], this.getGCPerc(gcCount, bpCount));
                                this.histograms.get(prevChr).writeToTemp();
                            }
                        }
                        skip = true;
                        prevChr = "NA";
                        continue;
                    }

                    this.histograms.put(line, new GCHistogram(line, rand));
                    if(prevChr.equals("NA")){
                        prevChr = line;
                        starts = wins.getStarts(line);
                        ends = wins.getEnds(line);
                    }else{
                        log.log(Level.FINEST, "Writing GC percentage of chr: " + prevChr);
                        if(binIdx < starts.length){
                            // Avoids issue where the binIdx is too high because of exact multiple of window size for the fasta entry
                            this.histograms.get(prevChr).addHistogram(prevChr, starts[binIdx], ends[binIdx], this.getGCPerc(gcCount, bpCount));
                            this.histograms.get(prevChr).writeToTemp();
                        }
                        prevChr = line;
                        totalCount = 0; bpCount = 0; gcCount = 0; binIdx = 0;
                        starts = wins.getStarts(line);
                        ends = wins.getEnds(line);
                    }
                }else if(!skip){
                    line = line.trim();
                    for(int i = 0; i < line.length(); i++){
                        bpCount++; totalCount++;
                        switch(line.charAt(i)){
                            case 'g':
                            case 'G':
                            case 'c':
                            case 'C':
                                gcCount++;
                        }
                        if(binIdx >= starts.length)
                            // Account for those smaller than readable regions
                            break;
                        if(bpCount >= wins.getWindowSize()){
                            this.histograms.get(prevChr).addHistogram(prevChr, starts[binIdx], ends[binIdx], this.getGCPerc(gcCount, bpCount));
                            bpCount = 0; gcCount = 0;
                            binIdx++;
                        }
                    }
                }                
            }
            // Add the last window
            if(binIdx < starts.length)
                this.histograms.get(prevChr).addHistogram(prevChr, starts[binIdx], ends[binIdx], this.getGCPerc(gcCount, bpCount));
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error reading fasta file: " + this.fastaPath.toString(), ex);
        }
        log.log(Level.INFO, "Created GC profiles for " + this.histograms.size() + " chromosomes.");
    }
    
    private double getGCPerc(int gcCount, int bpCount){
        if(bpCount == 0 || gcCount > bpCount)
            return 0.0d;
        return gcCount / (double) bpCount;
    }
    
    private boolean checkIfProfileFileExists(int windowsize){
        this.expectedFastaProfile = this.fastaPath.toString() +"." + windowsize + ".gc";
        return new File(this.expectedFastaProfile).canRead();
    }
    
    
    /*
        GC profile standard:
    TOP HEADER:
        1 byte  Magic byte
        4 bytes number of chrs
    HEADER:
        2 bytes chrname length
        X bytes chrname 
        4 bytes number of bins
    BINS:
        4 bytes start
        4 bytes end
        8 bytes gc percentage
    */
    private void generateGCHistoFromProfile(BamMetadataSampler bamMeta){
        // Load existing GC profile and transfer to histogram class
        log.log(Level.INFO, "Loading previously generated GC profile");
        ThreadTempRandAccessFile randprofile = new ThreadTempRandAccessFile(Paths.get(this.tmpPath.toString() + ".gcprofile.tmp"));
        try(RandomAccessFile rand = new RandomAccessFile(new File(this.expectedFastaProfile), "r")){
            // Get the top header
            byte magic = rand.readByte();
            byte[] ints = new byte[4];
            byte[] smInt = new byte[2];
            if(magic != 7)
                throw new Exception("Oops! The GC profile file might be corrupted!");
            
            rand.read(ints);
            
            int numchrs = IntUtils.byteArrayToInt(ints);
            log.log(Level.FINEST, "Predicted " + numchrs + " number of chromosomes in GC profile");
            for(int i = 0; i < numchrs; i++){
                rand.read(smInt);
                int chrnamelen = IntUtils.byteArrayToInt(smInt);
                
                byte[] chrs = new byte[chrnamelen];
                rand.read(chrs);
                String chrname = new String(chrs);
                
                if(!bamMeta.chrOrder.contains(chrname))
                    log.log(Level.WARNING, "Warning! GC profile contains chr not found in bam file: " + chrname);
                
                rand.read(ints);
                int numbins = IntUtils.byteArrayToInt(ints);
                this.histograms.put(chrname, new GCHistogram(chrname, randprofile));
                this.histograms.get(chrname).transferFromProfile(rand, numbins);
                log.log(Level.INFO, "Loaded GC profile from chr: " + chrname);
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error reading GC profile for fasta file: " + this.fastaPath.toString(), ex);
        }
        log.log(Level.INFO, "Loaded GC profile!");
    }
    
    private void writeGCProfileOut(BamMetadataSampler bamMeta){
        // Generate a GC profile index so that we don't have to do this again!
        log.log(Level.INFO, "Writing data to GC profile file: " + this.expectedFastaProfile);
        try(RandomAccessFile rand = new RandomAccessFile(new File(this.expectedFastaProfile), "rw")){
            // Writing out top header!
            rand.write(7);
            Long validChrs = this.histograms.keySet().stream().filter(s -> this.histograms.get(s).getNumEntries() > 0).count();
            rand.write(IntUtils.Int32ToByteArray(validChrs.intValue()));
            for(String chr : this.histograms.keySet()){
                GCHistogram gc = this.histograms.get(chr);
                gc.transferToProfile(rand);
                log.log(Level.INFO, "Finished writing to chr: " + chr);
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error writing out GC profile for fasta file: " + this.fastaPath.toString(), ex);
        }
    }
    
    public GCHistogram getGCHistogram(String chr){
        return this.histograms.get(chr);
    }
}
