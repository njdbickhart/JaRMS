/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils.ThreadingUtils;

import java.io.IOException;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Derek.Bickhart
 * @param <T>
 */
public class ThreadFactoryManager <T extends ThreadHistoFactory> {
    private static final Logger log = Logger.getLogger(ThreadFactoryManager.class.getName());
    
    private final int threads;
    private final ThreadTempRandAccessFile rand;
    
    public ThreadFactoryManager(int threads, ThreadTempRandAccessFile rand){
        this.threads = threads;
        this.rand = rand;
    }
        
    public void RunThreads(Map<T, Set<String>> grouping){
        ExecutorService executor = Executors.newFixedThreadPool(threads);
        
        for(Entry<T, Set<String>> entry : grouping.entrySet()){
            entry.getKey().ProcessSpecificWorkload(entry.getValue());
            executor.submit(entry.getKey());
            log.log(Level.FINE, "Submitted ThreadHistoFactory for processing: " + entry.getKey().getClass().getName());
        }
        
        executor.shutdown();
        
        while(!executor.isTerminated()){}
        
        log.log(Level.INFO, "Finished submission of factory pool. Consolidating temp files to: " + this.rand.GetFileName());
        for(T factory : grouping.keySet()){
            factory.Consolidate(rand);
        }
        try {
            rand.printIndex();
        } catch (IOException ex) {
            log.log(Level.SEVERE, "Error printing index file!", ex);
        }
    }
}
