/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jarms;

import GetCmdOpt.SimpleModeCmdLineParser;
import jarms.logger.ConsoleFormat;
import jarms.logger.LogFormat;
import jarms.modes.CallMode;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.ConsoleHandler;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author Derek.Bickhart
 */
public class JaRMS {
    private static final String version = "0.0.1";
    private static final Logger log = Logger.getLogger(JaRMS.class.getName());
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // Prepare and read command line options
        SimpleModeCmdLineParser cmd = PrepareCMDOptions();        
        cmd.GetAndCheckMode(args);
        
        // Set loggers and basic command options
        boolean debug = cmd.GetValue("debug").equals("true");
        setFileHandler(cmd.CurrentMode, args, debug);
        log.log(Level.INFO, "[MAIN] JaRMS version: " + version);
        String workingDir = System.getProperty("user.dir");
        log.log(Level.FINE, "[MAIN] Current working directory: " + workingDir);
        System.setProperty("java.io.tmpdir", workingDir);
        
        // Identify and run mode
        switch(cmd.CurrentMode){
            case "call":
                log.log(Level.INFO, "[MAIN] JaRMS call mode selected.");
                CallMode cluster = new CallMode(cmd);
                cluster.run();
                break;
            default:
                System.err.println("Error! Must designate a usage mode!");
                System.err.println(cmd.GetUsage());
                System.exit(-1);
        }
        
        log.log(Level.INFO, "[MAIN] Finished mode: " + cmd.CurrentMode);
        System.exit(0);
    }
    private static SimpleModeCmdLineParser PrepareCMDOptions(){
        String nl = System.lineSeparator();
        SimpleModeCmdLineParser cmd = new SimpleModeCmdLineParser("JaRMS\tA Java implementation of Mean-Shift RD signal decomposition" + nl
                + "Version: " + version + nl
            + "Usage: java -jar JaRMS.jar [mode] [mode specific options]" + nl
                + "Modes:" + nl
                + "\tcall\tCalls CNVs from BAM file" + nl,
                "call"
        );
        
        cmd.AddMode("call", 
                "JaRMS call mode" + nl +
                "Usage: java -jar JaRMS.jar cluster [-i bamfile -f fasta file -o output prefix] (option: -t number of threads)" + nl
                + "\t-i\tA BWA-processed bam file for processing" + nl
                + "\t-f\tThe reference genome fasta file that was used during alignment of the bam file" + nl
                + "\t-o\tOutput file prefix and directory" + nl
                + "\t-t\tNumber of threads to use [optional: use one thread]" + nl,
                "i:f:o:t:d|", 
                "ifo", 
                "ifotd", 
                "input", "fasta", "outbase", "threads", "debug");
        
        return cmd;
    }
    
    private static String loggerDate(){
        DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
        Date date = new Date();
        return dateFormat.format(date);
    }

    private static void setFileHandler(String type, String[] args, boolean debug) {
        // Create a log file and set levels for use with debugger
        FileHandler handler = null;
        ConsoleHandler console = null;
        String datestr = loggerDate();
        try {
            handler = new FileHandler("JaRMS." + type + "." + datestr + ".%u.%g.log");
            handler.setFormatter(new LogFormat());
            console = new ConsoleHandler();
            console.setFormatter(new ConsoleFormat());
            
            if(debug){
                handler.setLevel(Level.ALL);
                console.setLevel(Level.INFO);
            }else{
                handler.setLevel(Level.INFO);
                console.setLevel(Level.INFO);
            }
        } catch (IOException | SecurityException ex) {
            log.log(Level.SEVERE, "Error setting up logger!", ex);
        }
        Logger mainLog = Logger.getLogger("");
        // This will display all messages, but the handlers should filter the rest
        mainLog.setLevel(Level.ALL);
        for(Handler h : mainLog.getHandlers()){
            mainLog.removeHandler(h);
        }
        mainLog.addHandler(handler);
        mainLog.addHandler(console);
        
        // Log input arguments
        log.log(Level.INFO, "[MAIN] Command line arguments supplied: ");
        log.log(Level.INFO, StrUtils.StrArray.Join(args, " "));
        log.log(Level.INFO, "[MAIN] Debug flag set to: " + debug);
    }
    
}
