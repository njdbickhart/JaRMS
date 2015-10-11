/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author desktop
 */
public class ThreadTempRandAccessFile {
    private final Logger log = Logger.getLogger(ThreadTempRandAccessFile.class.getName());
    private final Map<String, Long> chrPointer = new ConcurrentHashMap<>();
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
    
    public void printIndex() throws IOException{
        try(BufferedWriter output = Files.newBufferedWriter(indexFile, Charset.defaultCharset())){
            this.chrPointer.entrySet().stream()
                    .forEach((s) -> {
                        try {
                            output.write(s.getKey() + "\t" + s.getValue() + System.lineSeparator());
                        } catch (IOException ex) {
                            log.log(Level.SEVERE, "Error writing to index file!", ex);
                        }
                    });
        }catch(IOException ex){
            
        }
    }
}
