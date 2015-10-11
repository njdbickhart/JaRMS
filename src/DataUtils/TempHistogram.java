/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Base temporary file class for storing and reloading histogram data into memory
 * @author Derek.Bickhart
 * @param <T>
 */
public abstract class TempHistogram <T extends Object>{
    protected ThreadTempRandAccessFile tempFile;
    protected int numEntries;
    protected String chr;
    protected int curIdx = 0;
    protected List<Integer> start = new ArrayList<>();
    protected List<Integer> end = new ArrayList<>();
    protected List<T> score = new ArrayList<>();
    
    public TempHistogram(String chr, ThreadTempRandAccessFile tmpFile){
        this.chr = chr;
        this.tempFile = tmpFile;
        //this.createTemp(tmpdir, chr);
    }
    
    protected void createTemp(Path path, String chr){
        /*try {
        Random rand = new Random();
        path = Paths.get(path.toString() + "." + chr + "." + rand.nextInt());
        //this.tempFile = Files.createTempFile(path.toString(), ".tmp");
        //this.tempFile.toFile().deleteOnExit();
        } catch (IOException ex) {
        Logger.getLogger(TempHistogram.class.getName()).log(Level.SEVERE, null, ex);
        }*/
    }
    
    
    
    public void resetIdx(){
        this.curIdx = 0;
    }
    
    /*
        Abstract methods
    */
    
    public abstract void addHistogram(String chr, int start, int end);
    
    public abstract void addHistogram(String chr, int start, int end, T score);
    
    /*
        Default byte block implementation:
        4 bytes: start
        4 bytes: end
        x bytes: score
    */    
    public abstract void writeToTemp();
    
    public abstract void readTemp();
    
    public abstract void clearData();
    
    public int iterateAdvance(){
        this.curIdx++;
        if(this.curIdx == this.numEntries){
            this.curIdx = 0;
            return -1;
        }
        return this.curIdx;
    }
    
    /*
        Getters
    */
    
    public int getCurIdx(){
        return this.curIdx;
    }
    
    public T getScore(){
        return this.score.get(curIdx);
    }
    
    public T getScore(int idx){
        return this.score.get(idx);
    }
    
    public int getEnd(){
        return this.end.get(curIdx);
    }
    
    public int getEnd(int idx){
        return this.end.get(idx);
    }
    
    public int getStart(){
        return this.start.get(curIdx);
    }
    
    public int getStart(int idx){
        return this.end.get(idx);
    }
    
    public String getChr(){
        return this.chr;
    }
    
    public int getNumEntries(){
        return this.numEntries;
    }
    
    /*
        Setters
    */
    
    public void setScore(T score){
        this.score.set(curIdx, score);
    }
    
    public void setScore(int idx, T score){
        this.score.set(idx, score);
    }
}
