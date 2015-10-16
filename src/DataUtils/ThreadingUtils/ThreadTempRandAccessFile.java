/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils.ThreadingUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author desktop
 */
public class ThreadTempRandAccessFile {
    private static final Logger log = Logger.getLogger(ThreadTempRandAccessFile.class.getName());
    private final Map<String, Long> chrPointer = new ConcurrentHashMap<>();
    private final Map<String, Long> chrLength = new ConcurrentHashMap<>();
    private RandomAccessFile dataFile;
    private final Path indexFile;
    
    public ThreadTempRandAccessFile(Path tempFile){
        try {
            this.dataFile = new RandomAccessFile(tempFile.toFile(), "rw");
        } catch (FileNotFoundException ex) {
            log.log(Level.SEVERE, "Could not create temp file!", ex);
        }
        this.indexFile = Paths.get(tempFile.toString() + ".idx");
    }
    
    public RandomAccessFile getFileForWriting(String chr){
        try {
            if(this.dataFile.getFilePointer() != this.dataFile.length()){
                this.dataFile.seek(this.dataFile.length());
            }
            this.chrPointer.put(chr, this.dataFile.length());
            
        } catch (IOException ex) {
            log.log(Level.SEVERE, "Error preparing file for writing!", ex);
        }
        
        return this.dataFile;
    }
    
    public RandomAccessFile getFileForReading(String chr){
        try {
            this.dataFile.seek(this.chrPointer.get(chr));
        } catch (IOException ex) {
            log.log(Level.SEVERE, "Error preparing file for writing!", ex);
        }
        return this.dataFile;
    }
    
    public void printIndex() {
        try(BufferedWriter output = Files.newBufferedWriter(indexFile, Charset.defaultCharset())){
            // Calculate chromosome long value length
            Set<Entry<String, Long>> pointers = this.chrPointer.entrySet();
            List<Entry<String, Long>> sortedSet = new ArrayList(pointers);
            Collections.sort(sortedSet, (e1, e2) -> e1.getValue().compareTo(e2.getValue()));
            
            String firstchr = sortedSet.get(0).getKey();
            long firstval = this.chrPointer.get(sortedSet.get(1).getKey());
            this.chrLength.put(firstchr, firstval);
            for(int i = 1; i < sortedSet.size() - 1; i++){
                this.chrLength.put(sortedSet.get(i).getKey(), this.chrPointer.get(sortedSet.get(i+1).getKey()) - sortedSet.get(i).getValue());
            }
            this.chrLength.put(sortedSet.get(sortedSet.size() - 1).getKey(), this.dataFile.length() - sortedSet.get(sortedSet.size() - 1).getValue());
            this.chrPointer.entrySet().stream()
                    .forEach((s) -> {
                        try {
                            output.write(s.getKey() + "\t" + s.getValue() + "\t" + this.chrLength.get(s.getKey()) + System.lineSeparator());
                        } catch (IOException ex) {
                            log.log(Level.SEVERE, "Error writing to index file!", ex);
                        }
                    });
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error writing to index file!" , ex);
        }
    }
    
    public boolean CanResume(){
        if(this.indexFile.toFile().canRead()){
            try {
                this.recoverFile();
            } catch (Exception ex) {
                log.log(Level.SEVERE, "Failed to resume file: " + this.GetFileName(), ex);
                return false;
            }
            return true;
        }else
            return false;
    }
    
    private void recoverFile(){
        try(BufferedReader input = Files.newBufferedReader(indexFile, Charset.defaultCharset())){
            String line;
            while((line = input.readLine()) != null){
                line = line.trim();
                String[] segs = line.split("\t");
                this.chrPointer.put(segs[0], Long.valueOf(segs[1]));
                this.chrLength.put(segs[0], Long.valueOf(segs[2]));
            }
        }catch(IOException ex){
            log.log(Level.SEVERE, "Error recovering temp file from index: " + this.indexFile.toString(), ex);
        }
    }
    
    public void Close(){
        try {
            this.dataFile.close();
            this.indexFile.toFile().delete();
        } catch (IOException ex) {
            log.log(Level.SEVERE, "Could not close temp file!", ex);
        }
        
    }
    
    public String GetFileName(){
        return this.indexFile.toString();
    }
    
    public Set<String> getListChrs(){
        return this.chrPointer.keySet();
    }
    
    public long getChrLength(String chr){
        if(this.chrLength.containsKey(chr))
            return this.chrLength.get(chr);
        else return 0l;
    }
}
