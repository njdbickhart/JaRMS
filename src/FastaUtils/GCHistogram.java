/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaUtils;

import DataUtils.ThreadingUtils.TempHistogram;
import DataUtils.ThreadingUtils.ThreadTempRandAccessFile;
import TempFiles.binaryUtils.DoubleUtils;
import TempFiles.binaryUtils.IntUtils;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.nio.file.Path;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Derek.Bickhart
 */
public class GCHistogram extends TempHistogram<Double>{
    private static final Logger log = Logger.getLogger(GCHistogram.class.getName());
    
    public GCHistogram(String chr, ThreadTempRandAccessFile tmpdir) {
        super(chr, tmpdir);
    }

    @Override
    public void addHistogram(String chr, int start, int end) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void addHistogram(String chr, int start, int end, Double score) {
        super.start.add(start);
        super.end.add(end);
        super.score.add(score);
        super.numEntries++;
    }

    @Override
    public void writeToTemp() {
        try{
            RandomAccessFile rand = this.tempFile.getFileForWriting(chr);
            for(int x = 0; x < super.numEntries; x++){
                rand.write(IntUtils.Int32ToByteArray(this.start.get(x)));
                rand.write(IntUtils.Int32ToByteArray(this.end.get(x)));
                rand.write(DoubleUtils.Dbl64toByteArray(this.score.get(x)));
            }
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error writing to GCHistogram temp file for chr: " + this.chr, ex);
        }
        
        // Clear the list
        this.clearData();
    }
    
    public void transferFromProfile(RandomAccessFile profile, int numBins){
        try{
            RandomAccessFile rand = this.tempFile.getFileForWriting(chr);
            byte[] ints = new byte[4];
            byte[] dbls = new byte[8];
            this.numEntries = numBins;
            
            for(int i = 0; i < numBins; i++){
                profile.read(ints);
                rand.write(ints);
                
                profile.read(ints);
                rand.write(ints);
                
                profile.read(dbls);
                rand.write(dbls);
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error transferring profile for chr: " + this.chr, ex);
        }
    }
    
    /*
        GC profile standard:
    HEADER:
        2 bytes chrname length
        X bytes chrname 
        4 bytes number of bins
    BINS:
        4 bytes start
        4 bytes end
        8 bytes gc percentage
    */
    public void transferToProfile(RandomAccessFile profile) throws Exception{
        try{
            RandomAccessFile rand = this.tempFile.getFileForReading(chr);
            // Start with writing the header
            // Account for super long contig names!
            byte[] chrlen = IntUtils.Int16ToTwoByteArray(this.chr.length());
            byte[] chrname = this.chr.getBytes();
            byte[] numbins = IntUtils.Int32ToByteArray(this.numEntries);
            profile.write(chrlen);
            profile.write(chrname);
            profile.write(numbins);
            
            // Now, start writing the bins out
            byte[] ints = new byte[4];
            byte[] dbls = new byte[8];
            
            for(int i = 0; i < this.numEntries; i++){
                rand.read(ints);
                profile.write(ints);
                
                rand.read(ints);
                profile.write(ints);
                
                rand.read(dbls);
                profile.write(dbls);
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error transferring profile for chr: " + this.chr, ex);
        }
    }

    @Override
    public void readTemp() {
        try{
            RandomAccessFile rand = this.tempFile.getFileForReading(chr);
            if(this.numEntries <= 0)
                throw new Exception ("Reading empty temp file!");
            byte[] ints = new byte[4];
            byte[] dbls = new byte[8];
            for(int x = 0; x < super.numEntries; x++){
                rand.read(ints);
                this.start.add(IntUtils.byteArrayToInt(ints));
                
                rand.read(ints);
                this.end.add(IntUtils.byteArrayToInt(ints));
                
                rand.read(dbls);
                this.score.add(DoubleUtils.BytetoDouble(dbls));
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error reading from GCHistogram temp file for chr: " + this.chr, ex);
        }
    }
    
    // GC score is going to be from 0 - 1000 to account for array index values
    public int getRoundedGC(int idx){
        return (int)(this.score.get(idx) * 1000.0d);
    }
    public int getRoundedGC(){
        return (int)(this.score.get(this.curIdx) * 1000.0d);
    }
    
    public double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }
    
    @Override
    public void clearData() {
        this.start.clear();
        this.end.clear();
        this.score.clear();
        
        System.gc();
    }
}
