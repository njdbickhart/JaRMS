/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils.ThreadingUtils;

import java.util.Set;

/**
 *
 * @author Derek.Bickhart
 */
public interface ThreadHistoFactory extends Runnable{
    public void ProcessSpecificWorkload(Set<String> chrs);
    public void Consolidate(ThreadTempRandAccessFile rand);
    public void ResumeFromTempFile(ThreadTempRandAccessFile rand);
}
