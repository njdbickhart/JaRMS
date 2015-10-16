/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

import DataUtils.ThreadingUtils.TempHistogram;
import DataUtils.ThreadingUtils.ThreadHistogram;
import DataUtils.ThreadingUtils.ThreadTempRandAccessFile;
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
public class ChrHistogram extends TempHistogram<Double> implements ThreadHistogram{
    private static final Logger log = Logger.getLogger(ChrHistogram.class.getName());
    private double sum = 0.0d;
    
    public ChrHistogram(String chr, ThreadTempRandAccessFile tmpdir) {
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
        this.sum += score;
        if(this.sum == Double.NaN){
            System.err.println("Error with sum score!");
        }
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
        try{
            RandomAccessFile rand = this.tempFile.getFileForWriting(chr);
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
            log.log(Level.SEVERE, "Error reading from ChrHistogram temp file for chr: " + this.chr, ex);
        }
    }
    
    public double[] retrieveRDBins(){
        double[] bins = new double[this.numEntries];
        try{
            RandomAccessFile rand = this.tempFile.getFileForReading(chr);
            if(this.numEntries <= 0)
                throw new Exception ("Reading empty temp file!");
            byte[] ints = new byte[4];
            byte[] dbls = new byte[8];
            for(int x = 0; x < super.numEntries; x++){
                rand.read(ints);
                
                rand.read(ints);
                
                rand.read(dbls);
                bins[x] = DoubleUtils.BytetoDouble(dbls);
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error reading bins from ChrHistogram temp file for chr: " + this.chr, ex);
        }
        
        return bins;
    }

    @Override
    public void clearData() {
        this.start.clear();
        this.end.clear();
        this.score.clear();
        
        System.gc();
    }
    
    
    public double getMean(){
        return this.sum / this.numEntries;
    }
    
    public double getSum(){
        log.log(Level.FINEST, "Sum score for chr: " + this.chr + " " + this.sum);
        return this.sum;
    }
    
    public double getSumSquares(double mean){
        double ss = 0.0d;
        try{
            RandomAccessFile rand = this.tempFile.getFileForReading(chr);
            if(this.numEntries <= 0)
                throw new Exception ("Reading empty temp file!");
            byte[] ints = new byte[4];
            byte[] dbls = new byte[8];
            for(int x = 0; x < super.numEntries; x++){
                rand.read(ints);
                
                rand.read(ints);
                
                rand.read(dbls);
                double value = DoubleUtils.BytetoDouble(dbls);
                
                ss += Math.pow(value - mean, 2.0d);
            }
        }catch(Exception ex){
            log.log(Level.SEVERE, "Error reading bins from ChrHistogram temp file for chr: " + this.chr, ex);
        }
        return ss;
    }
    
    public void recalculateSumScore(){
        if(Double.isNaN(sum)){
            try{
                RandomAccessFile rand = this.tempFile.getFileForReading(chr);
                if(this.numEntries <= 0)
                    throw new Exception ("Reading empty temp file!");
                byte[] ints = new byte[4];
                byte[] dbls = new byte[8];
                for(int x = 0; x < super.numEntries; x++){
                    rand.read(ints);

                    rand.read(ints);

                    rand.read(dbls);
                    double value = DoubleUtils.BytetoDouble(dbls);
                    if(Double.isNaN(value)){
                        value = 0.0d;
                    }
                    this.sum += value;
                }
            }catch(Exception ex){
                log.log(Level.SEVERE, "Error reading bins from ChrHistogram temp file for chr: " + this.chr, ex);
            }
        }
    }

    @Override
    public void TransferToSharedTemp(ThreadTempRandAccessFile rand) {
        try{
            byte[] ints = new byte[4];
            byte[] dbls = new byte[8];
            RandomAccessFile thisTemp = this.tempFile.getFileForReading(chr);
            RandomAccessFile temp = rand.getFileForWriting(chr);
            for(int x = 0; x < super.numEntries; x++){
                temp.write(thisTemp.read(ints));
                temp.write(thisTemp.read(ints));
                temp.write(thisTemp.read(dbls));
            }
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error transfering ChrHistogram temp file for chr: " + this.chr, ex);
        }
        
    }
    
    public void setNumEntries(int value){
        this.numEntries = value;
    }
}
