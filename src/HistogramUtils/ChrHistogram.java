/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

import DataUtils.TempHistogram;
import TempFiles.binaryUtils.DoubleUtils;
import TempFiles.binaryUtils.IntUtils;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Path;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Derek.Bickhart
 */
public class ChrHistogram extends TempHistogram<Double>{
    private static final Logger log = Logger.getLogger(ChrHistogram.class.getName());
    private boolean[] freeze;
    
    public ChrHistogram(String chr, Path tmpdir) {
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

    /*
        Default byte block implementation:
        4 bytes: start
        4 bytes: end
        8 bytes: score
    */
    @Override
    public void writeToTemp() {
        try(RandomAccessFile rand = new RandomAccessFile(this.tempFile.toFile(), "rw")){
            for(int x = 0; x < super.numEntries; x++){
                rand.write(IntUtils.Int32ToByteArray(this.start.get(x)));
                rand.write(IntUtils.Int32ToByteArray(this.end.get(x)));
                rand.write(DoubleUtils.Dbl64toByteArray(this.score.get(x)));
            }
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error writing to ChrHistogram temp file for chr: " + this.chr, ex);
        }
        
        // Clear the list
        this.clearData();
    }

    @Override
    public void readTemp() {
        try(RandomAccessFile rand = new RandomAccessFile(this.tempFile.toFile(), "r")){
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
            log.log(Level.SEVERE, "Error reading from ChrHistogram temp file for chr: " + this.chr, ex);
        }
    }

    @Override
    public void clearData() {
        this.start.clear();
        this.end.clear();
        this.score.clear();
        
        System.gc();
    }
    
}
