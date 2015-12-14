/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jarms.modes;

import DataUtils.ThreadingUtils.ThreadDivisor;
import DataUtils.ThreadingUtils.ThreadFactoryManager;
import DataUtils.ThreadingUtils.ThreadTempRandAccessFile;
import DataUtils.WindowPlan;
import FastaUtils.GlobalGCCorrectionProfile;
import FastaUtils.HTSGCWindowFactory;
import GetCmdOpt.ArrayModeCmdLineParser;
import HistogramUtils.BamMetadataSampler;
import HistogramUtils.ChrHistogramFactory;
import SVCaller.MeanShiftMethod;
import SVCaller.SVSegmentation;
import htsjdk.samtools.SamReader;
import java.io.File;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author desktop
 */
public class CallMode {
    private static final Logger log = Logger.getLogger(CallMode.class.getName());
    private final ArrayModeCmdLineParser cmd;
    
    private final List<String> bamFiles;
    private final String outDir;
    private final String fastaFile;
    private int OverrideWinSize = -1;
    private int threads = 1;

    public CallMode(ArrayModeCmdLineParser cmd) {
        this.cmd = cmd;
        
        // Checking to see if we can read at least one Bam file!
        if(cmd.HasOpt("input") && new File(cmd.GetArray("input").get(0)).canRead())
            this.bamFiles = cmd.GetArray("input");
        else{
            this.bamFiles = null;
            ErrorExit("input");
        }
        
        if(cmd.HasOpt("outbase"))
            this.outDir = cmd.GetValue("outbase");
        else
            this.outDir = null;
        
        if(cmd.HasOpt("fasta") && new File(cmd.GetValue("fasta")).canRead())
            this.fastaFile = cmd.GetValue("fasta");
        else{
            this.fastaFile = null;
            ErrorExit(cmd.GetValue("fasta"));
        }
        
        if(cmd.HasOpt("window"))
            this.OverrideWinSize = Integer.parseInt(cmd.GetValue("window"));
        
        if(cmd.HasOpt("threads"))
            this.threads = Integer.parseInt(cmd.GetValue("threads"));
    }

    public void run() {
        // Get BAM metadata
        BamMetadataSampler metadata = new BamMetadataSampler(this.bamFiles);
        metadata.getMetaData();
        
        // Identify window plan
        WindowPlan wins = new WindowPlan();
        wins.GenerateWindows(metadata, this.OverrideWinSize);
        
        //if(this.threads == 1)
            SingleThreadRun(metadata, wins);
        //else if(this.threads > 1)
            //MultiThreadRun(metadata, wins);
    }
    
    private void MultiThreadRun(BamMetadataSampler metadata, WindowPlan wins){
        // Divide tasks
        List<Set<String>> chrs = ThreadDivisor.EstimateThreadDivision(threads, metadata.chrOrder);
        
        // Create RD histograms
        log.log(Level.FINE, "[CALLMODE] Starting RD histogram counting, multi-threaded");
        Map<ChrHistogramFactory, Set<String>> rdworkers = new HashMap<>();
        for(int i = 0; i < chrs.size(); i++){
            SamReader reader = null;
            /*SamReaderFactory.make()
            .validationStringency(ValidationStringency.LENIENT)
            .samRecordFactory(DefaultSAMRecordFactory.getInstance())
            .open(new File());*/
            ThreadTempRandAccessFile threadtemp = new ThreadTempRandAccessFile(Paths.get(outDir + "." + i + ".rdhisto.tmp"));
            rdworkers.put(new ChrHistogramFactory(reader, wins, threadtemp), chrs.get(i));
        }
        ThreadTempRandAccessFile rdtemp = new ThreadTempRandAccessFile(Paths.get(outDir + ".rdhisto.tmp"));
        ThreadFactoryManager rdManager = new ThreadFactoryManager(this.threads, rdtemp);
        rdManager.RunThreads(rdworkers);
        
        log.log(Level.FINE, "[CALLMODE] Finished RD histogram start");
        
    }

    private void SingleThreadRun(BamMetadataSampler metadata, WindowPlan wins) {
        
        // Create RD histograms
        log.log(Level.FINE, "[CALLMODE] Starting RD histogram counting");
        ThreadTempRandAccessFile rawHistoRand = new ThreadTempRandAccessFile(Paths.get(this.outDir + ".rdhisto.tmp"));
        ChrHistogramFactory rawRDHisto = new ChrHistogramFactory(rawHistoRand);
        if(!rawHistoRand.CanResume()){
            try {
                rawRDHisto.processMultipleBamNoRG(bamFiles, wins, outDir);
                //rawRDHisto.processBamNoRG(bamFile, wins, outDir);
            } catch (Exception ex) {
                log.log(Level.SEVERE, "[CALLMODE] Error processing bam file!", ex);
            }
            log.log(Level.FINE, "[CALLMODE] Finished RD histogram start");
            rawHistoRand.printIndex();
        }else{
            log.log(Level.FINE, "[CALLMODE] Resumed RD calculation from previous temp files");
            rawRDHisto.ResumeFromTempFile(rawHistoRand);
        }
        // Calculate global mean and SD RD prior to correction
        double rawSum = wins.getChrStream()
                .filter(s -> rawRDHisto.hasChrHistogram(s))
                .map(s -> rawRDHisto.getChrHistogram(s).getSum())
                .reduce(0.0d, (a, b) -> a + b);
        int rawNumBins = wins.getChrStream()
                .filter(s -> rawRDHisto.hasChrHistogram(s))
                .map(s -> rawRDHisto.getChrHistogram(s).getNumEntries())
                .reduce(0, (a, b) -> a + b);
        final double rawMean = rawSum / (double) rawNumBins;
        
        double rawSS = wins.getChrStream()
                .filter(s -> rawRDHisto.hasChrHistogram(s))
                .map(s -> rawRDHisto.getChrHistogram(s).getSumSquares(rawMean))
                .reduce(0.0d, (a, b) -> a + b);
        rawSS /= rawNumBins - 1;
        
        final double rawStdev = Math.sqrt(rawSS);
        log.log(Level.INFO, "[CALLMODE] Unadjusted RD values: mean-> " + rawMean + " sd-> " +rawStdev);
        
        ThreadTempRandAccessFile gcCorrRand = new ThreadTempRandAccessFile(Paths.get(this.outDir + ".gccorr.tmp"));
        if(new File(gcCorrRand.GetFileName()).exists()){
            log.log(Level.INFO, "[CALLMODE] Found index file: " + gcCorrRand.GetFileName());
        }else{
            log.log(Level.INFO, "[CALLMODE] Could not find index file: " + gcCorrRand.GetFileName());
        }
        
        ChrHistogramFactory gcCorrectRDHisto = new ChrHistogramFactory(gcCorrRand);
        if(gcCorrRand.CanResume()){
            log.log(Level.FINE, "[CALLMODE] Resumed GC calculation from previous temp files");
            gcCorrectRDHisto.ResumeFromTempFile(gcCorrRand);
        }else{
            // Generate GC correction scheme
            log.log(Level.FINE, "[CALLMODE] Calculating GC windows");
            ThreadTempRandAccessFile gcProfilerand = new ThreadTempRandAccessFile(Paths.get(this.outDir + ".gcprofile.tmp"));
            HTSGCWindowFactory GCWins = new HTSGCWindowFactory(this.fastaFile, gcProfilerand);
            GCWins.generateGCProfile(metadata, wins);
            log.log(Level.FINE, "[CALLMODE] Estimated GC profile");

            // Create GC correction utility

            GlobalGCCorrectionProfile gccorrect = new GlobalGCCorrectionProfile();
            gccorrect.CalculationGCCorrectionValues(metadata, rawRDHisto, GCWins, rawMean);
            log.log(Level.FINE, "[CALLMODE] Corrected GC bias");

            // Correct GC
            gcCorrectRDHisto = gccorrect.CorrectGC(gcCorrRand, metadata, rawRDHisto, GCWins);

            // Check to see if we nee to recalculate the sum values for each GCCorrected entity
            gcCorrectRDHisto.checkSumScores();
            gcCorrRand.printIndex();
        }
        
        
        // Mean shift signal
        log.log(Level.FINE, "[CALLMODE] Beginning mean shift algorithm");
        ThreadTempRandAccessFile levelRand = new ThreadTempRandAccessFile(Paths.get(this.outDir + ".levels.tmp"));
        MeanShiftMethod shifter = new MeanShiftMethod();
        shifter.Partition(gcCorrectRDHisto, wins, levelRand, 128, threads);
        log.log(Level.FINE, "[CALLMODE] Ended mean shift algorithm");
        
        // Call SVs
        SVSegmentation svCaller = new SVSegmentation(this.outDir);
        svCaller.RunSegmentation(gcCorrectRDHisto, wins, shifter, threads);
        log.log(Level.FINE, "[CALLMODE] Called SVs");
        
        // Print out the calls and levels
        svCaller.printOutAllCalls();
        svCaller.printOutCondensedLevels(shifter, wins, svCaller.getGlobalMean());
        
        // Clean up tmp files
        final File folder = new File(System.getProperty("user.dir"));
        final File[] files = folder.listFiles((final File dir, final String name) -> 
                name.endsWith("gccorr.tmp") || name.endsWith("gcprofile.tmp") || name.endsWith("levels.tmp") || name.endsWith("rdhisto.tmp")
                || name.endsWith("rdhisto.tmp.idx") || name.endsWith("gcprofile.tmp.idx"));
        for ( final File file : files ) {
            if ( !file.delete() ) {
                System.err.println( "Can't remove " + file.getAbsolutePath() );
            }
        }
    }
    
    private void ErrorExit(String value){
        log.log(Level.SEVERE, "[CALLMODE] Error with option: " + value + " Could not find the file/directory!");
        System.exit(-1);
    }
}
