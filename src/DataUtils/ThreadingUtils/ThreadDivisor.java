/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils.ThreadingUtils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author Derek.Bickhart
 */
public class ThreadDivisor {
    public static List<Set<String>> EstimateThreadDivision(int threads, Set<String> chrs){
        int ratio = (int) Math.ceil(chrs.size() / (double)threads);
        List<Set<String>> tasks = new ArrayList<>(ratio);
        
        // Split up initial threads
        Object[] chrarray = chrs.toArray();
        int count = 0;
        for(int x = 0; x < chrs.size() - ratio; x += ratio){
            Set<String> temp = new HashSet<>();
            for(int y = x; y < x + ratio; y++){
                temp.add((String) chrarray[y]);
                count++;
            }
            tasks.add(temp);
        }
        
        // Catch remaining and append to final thread division
        if(count < chrs.size()){
            for(int y = count; y < chrs.size(); y++){
                tasks.get(tasks.size() - 1).add((String) chrarray[y]);
            }
        }
        
        return tasks;
    }
}
