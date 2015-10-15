/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaUtils;

import DataUtils.ThreadTempRandAccessFile;
import DataUtils.WindowPlan;
import HistogramUtils.BamMetadataSampler;
import TempFiles.binaryUtils.IntUtils;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import java.io.File;
import java.io.RandomAccessFile;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Derek.Bickhart
 */
public class HTSGCWindowFactory extends GCWindowFactory{
    private static final Logger log = Logger.getLogger(HTSGCWindowFactory.class.getName());
    
    //private final Map<String, GCHistogram> histograms = new HashMap<>();
    //private final Path fastaPath;
    //private final Path tmpPath;
    //private String expectedFastaProfile;

    public HTSGCWindowFactory(String fastaFile, String tmpdir){
        super(fastaFile, tmpdir);
        //this.fastaPath = Paths.get(fastaFile);
        //this.tmpPath = Paths.get(tmpdir);
    }
    
    @Override
    public void generateGCProfile(BamMetadataSampler bamMeta, WindowPlan wins){
        // Check if we can load existing profile file
        if(this.checkIfProfileFileExists(wins.getWindowSize())){
            this.generateGCHistoFromProfile(bamMeta);
        }else{
            this.processFastaFile(wins);
            this.writeGCProfileOut(bamMeta);
        }        
    }
    
    private boolean checkIfProfileFileExists(int windowsize){
        this.expectedFastaProfile = this.fastaPath.toString() +"." + windowsize + ".gc";
        return new File(this.expectedFastaProfile).canRead();
    }
    
    private void processFastaFile(WindowPlan wins){
        ThreadTempRandAccessFile rand = new ThreadTempRandAccessFile(Paths.get(this.tmpPath.toString() + ".gcprofile.tmp"));
        ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(this.fastaPath.toFile());
        
        for(String chr : wins.getChrList()){
            Integer[] starts = wins.getStarts(chr);
            Integer[] ends = wins.getEnds(chr);
            this.histograms.put(chr, new GCHistogram(chr, rand));
            
            log.log(Level.FINEST, "Calculating GC percentage for chr: " + chr);
            ReferenceSequence refSeq = refFile.getSequence(chr);
            byte[] seq = refSeq.getBases();
            
            for(int x = 0; x < starts.length; x++){
                int gccount = 0, bpcount = 0, ncount = 0;
                for(int y = starts[x]; y < ends[x] - 1 && y < seq.length; y++){
                    switch(seq[y]){
                        case 'G':
                        case 'C':
                        case 'c':
                        case 'g':
                            gccount++;
                            break;
                        case 'N':
                        case 'n':
                            ncount++;
                            break;
                    }
                }
                bpcount = ends[x] - starts[x];
                double gc = this.getGCPerc(gccount, bpcount);
                this.histograms.get(chr).addHistogram(chr, starts[x], ends[x], gc);
            }
            this.histograms.get(chr).writeToTemp();
        }
    }
    
    private double getGCPerc(int gcCount, int bpCount){
        if(bpCount == 0 || gcCount > bpCount)
            return 0.0d;
        return gcCount / (double) bpCount;
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
                log.log(Level.FINEST, "Loaded GC profile from chr: " + chrname);
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
    
    @Override
    public GCHistogram getGCHistogram(String chr){
        return this.histograms.get(chr);
    }
}
