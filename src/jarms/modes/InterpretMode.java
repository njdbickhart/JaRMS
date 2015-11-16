/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jarms.modes;

import DataUtils.ThreadingUtils.ThreadTempRandAccessFile;
import GetCmdOpt.ArrayModeCmdLineParser;
import HistogramUtils.ChrHistogramFactory;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Derek.Bickhart
 */
public class InterpretMode {
    private static final Logger log = Logger.getLogger(InterpretMode.class.getName());
    
    private final ArrayModeCmdLineParser cmd;
    private File binary;
    private File index;
    private Path output;
    
    public InterpretMode(ArrayModeCmdLineParser cmd){
        this.cmd = cmd;
        
        if(cmd.HasOpt("input")){
            this.binary = new File(cmd.GetValue("input"));
            this.index = new File(cmd.GetValue("input") + ".idx");
            
            if(!this.binary.canRead() || !this.index.canRead()){
                log.log(Level.SEVERE, "Error! Cannot read index or binary input files!");
                System.exit(-1);
            }
        }
        
        if(cmd.HasOpt("output")){
            this.output = Paths.get(cmd.GetValue("output"));
        }
    }
    
    public void run(){
        log.log(Level.FINE, "[INTERPRET] Starting RD histogram resumption");
        ThreadTempRandAccessFile HistoRand = new ThreadTempRandAccessFile(Paths.get(this.binary.getPath()));
        ChrHistogramFactory RDHisto = new ChrHistogramFactory(HistoRand);
        if(HistoRand.CanResume()){
            RDHisto.ResumeFromTempFile(HistoRand);
            try (BufferedWriter outfile = Files.newBufferedWriter(output, Charset.defaultCharset())){
                RDHisto.PrintWindowsFromTempFile(outfile);
            } catch (Exception ex) {
                log.log(Level.SEVERE, "[INTERPRET] Error processing binary file!", ex);
            }
            log.log(Level.FINE, "[INTERPRET] Finished RD histogram start");
        }
        log.log(Level.INFO, "[INTERPRET] Finished printing data to text");
    }
}
