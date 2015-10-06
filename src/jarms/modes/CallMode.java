/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jarms.modes;

import DataUtils.WindowPlan;
import FastaUtils.GCWindowFactory;
import FastaUtils.GlobalGCCorrectionProfile;
import GetCmdOpt.SimpleModeCmdLineParser;
import HistogramUtils.BamMetadataSampler;
import HistogramUtils.ChrHistogramFactory;
import SVCaller.MeanShiftMethod;
import SVCaller.SVSegmentation;
import java.io.File;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author desktop
 */
public class CallMode {
    private static final Logger log = Logger.getLogger(CallMode.class.getName());
    private final SimpleModeCmdLineParser cmd;
    
    private final String bamFile;
    private final String outDir;
    private final String fastaFile;
    private int threads = 1;

    public CallMode(SimpleModeCmdLineParser cmd) {
        this.cmd = cmd;
        
        if(cmd.HasOpt("input") && new File(cmd.GetValue("input")).canRead())
            this.bamFile = cmd.GetValue("input");
        else{
            this.bamFile = null;
            ErrorExit(cmd.GetValue("input"));
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
        
        if(cmd.HasOpt("threads"))
            this.threads = Integer.parseInt(cmd.GetValue("threads"));
    }

    public void run() {
        // Get BAM metadata
        BamMetadataSampler metadata = new BamMetadataSampler(this.bamFile);
        metadata.getMetaData();
        
        // Identify window plan
        WindowPlan wins = new WindowPlan();
        wins.GenerateWindows(metadata);
        
        // Create RD histograms
        ChrHistogramFactory rawRDHisto = new ChrHistogramFactory();
        try {
            rawRDHisto.processBamNoRG(bamFile, wins, outDir);
        } catch (Exception ex) {
            log.log(Level.SEVERE, "[CALLMODE] Error processing bam file!", ex);
        }
        // Calculate global mean and SD RD prior to correction
        double rawSum = wins.getChrList().stream()
                .map(s -> rawRDHisto.getChrHistogram(s).getSum())
                .reduce(0.0d, (a, b) -> a + b);
        int rawNumBins = wins.getChrList().stream()
                .map(s -> rawRDHisto.getChrHistogram(s).getNumEntries())
                .reduce(0, (a, b) -> a + b);
        final double rawMean = rawSum / (double) rawNumBins;
        
        double rawSS = wins.getChrList().stream()
                .map(s -> rawRDHisto.getChrHistogram(s).getSumSquares(rawMean))
                .reduce(0.0d, (a, b) -> a + b);
        rawSS /= rawNumBins - 1;
        
        final double rawStdev = Math.sqrt(rawSS);
        log.log(Level.INFO, "[CALLMODE] Unadjusted RD values: mean-> " + rawMean + " sd-> " +rawStdev);
        
        // Generate GC correction scheme
        GCWindowFactory GCWins = new GCWindowFactory(this.fastaFile, this.outDir);
        GCWins.generateGCProfile(metadata, wins);
        
        // Create GC correction utility
        GlobalGCCorrectionProfile gccorrect = new GlobalGCCorrectionProfile();
        gccorrect.CalculationGCCorrectionValues(metadata, rawRDHisto, GCWins, rawMean);
        
        // Correct GC
        ChrHistogramFactory gcCorrectRDHisto = gccorrect.CorrectGC(Paths.get(this.outDir), metadata, rawRDHisto, GCWins);
        
        // Mean shift signal
        MeanShiftMethod shifter = new MeanShiftMethod();
        shifter.Partition(rawRDHisto, wins, Paths.get(this.outDir), 128, threads);
        
        // Call SVs
        SVSegmentation svCaller = new SVSegmentation(this.outDir);
        svCaller.RunSegmentation(gcCorrectRDHisto, wins, shifter, threads);
        
        // Print out the calls
        svCaller.printOutAllCalls();
    }
    
    private void ErrorExit(String value){
        log.log(Level.SEVERE, "[CALLMODE] Error with option: " + value + " Could not find the file/directory!");
        System.exit(-1);
    }
}
