/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package HistogramUtils;

import DataUtils.ThreadTempRandAccessFile;
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
public class LevelHistogram extends ChrHistogram{
    private final Logger log = Logger.getLogger(LevelHistogram.class.getName());
    private final double PRECISION = 0.01d;

    public LevelHistogram(String chr, ThreadTempRandAccessFile tmpdir) {
        super(chr, tmpdir);
    }
    
    /*
        Method designed to rewrite previously written Level data and replace with the merged data
    */
    public void replaceTemp(){
        try{
            RandomAccessFile rand = this.tempFile.getFileForReading(chr);
            for(int x = 0; x < super.numEntries; x++){
                //rand.write(IntUtils.Int32ToByteArray(this.start.get(x)));
                //rand.write(IntUtils.Int32ToByteArray(this.end.get(x)));
                rand.seek(rand.getFilePointer() + 8);
                rand.write(DoubleUtils.Dbl64toByteArray(this.score.get(x)));
            }
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error writing to LevelHistogram temp file for chr: " + this.chr, ex);
        }
        
        // Clear the list
        this.clearData();
    }
    
    public void performMerger(double delta){
        if(this.score.isEmpty())
            this.readTemp();
        this.mergeLevels(0, delta);
    }
    
    private void mergeLevels(int recursion, double delta){
        if(recursion >= 1000000){
            return;
            // Just to prevent cases where we cannot smooth the signal beyond a large number of retries!
        }
        boolean ret = false;

        int start1 = 0, end1 = 0, start2, end2;
        while (end1 < this.numEntries &&
               Math.abs(this.score.get(end1) - this.score.get(start1)) < PRECISION) 
            end1++;
        end1--; 
        start2 = end2 = end1 + 1;
        while (end2 < this.numEntries &&
               Math.abs(this.score.get(end2) - this.score.get(start2)) < PRECISION) 
            end2++;
        end2--;

        while (start2 < this.numEntries) {
            double v1 = Math.abs(this.score.get(start1) - this.score.get(start2));
            if (v1 < delta) {
                double v2 = v1 + 1, v3 = v1 + 1;
                if (start1 - 1 >= 0)
                  v2 = Math.abs(this.score.get(start1) - this.score.get(start1 - 1));
                if (end2 + 1 < this.numEntries)
                  v3 = Math.abs(this.score.get(end2)   - this.score.get(end2 + 1));
                if (v1 < v2 && v1 < v3) {
                    ret = true;
                    double nl = (end1 - start1 + 1) * this.score.get(start1);
                    nl       += (end2 - start2 + 1) * this.score.get(start2);
                    nl       *= 1.0d / (end2 - start1 + 1);
                    for (int i = start1;i <= end2;i++) 
                        this.score.set(i, nl);
                    start2 = start1;
                }
            }
            start1 = start2;
            end1   = end2;
            start2 = end2 = end1 + 1;
            while (end2 < this.numEntries &&
                   Math.abs(this.score.get(end2) - this.score.get(start2)) < PRECISION) 
                end2++;
            end2--;
        }

        if(ret)
            this.mergeLevels(++recursion, delta);
    }
}
