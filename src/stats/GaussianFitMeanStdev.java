/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import DataUtils.WindowPlan;
import HistogramUtils.ChrHistogramFactory;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

/**
 *
 * @author Derek.Bickhart
 */
public class GaussianFitMeanStdev {
    private static Logger log = Logger.getLogger(GaussianFitMeanStdev.class.getName());
    
    GaussianCurveFitter fitter;
    private double fmean = 0.0d;
    private double fstdev = 0.0d;
    private double mean;
    private double stdev;
    
    public boolean cantUse = false;
    
    public void CalculateGlobalMeanStdev(ChrHistogramFactory gchisto, WindowPlan wins){
        Long maxvalue = 0l;
        int count = 0;
        for(String chr : wins.getChrList()){
            if(!gchisto.hasChrHistogram(chr))
                continue;
            for(Double d : gchisto.getChrHistogram(chr).retrieveRDBins()){
                count++;
                if(d.longValue() > maxvalue)
                    maxvalue = d.longValue();
            }
        }
        if(count < 3 || maxvalue.intValue() < 3){
            log.log(Level.SEVERE, "Error! Could not calculate Global Mean and Stdev! Count: " + count + " maxvalue: " + maxvalue);
            System.exit(-1);
        }
        
        Double[] bins = new Double[maxvalue.intValue() + 1];
        java.util.Arrays.fill(bins, 0.0d);
        double sum = 0.0d;
        for(String chr : wins.getChrList()){
            if(!gchisto.hasChrHistogram(chr))
                continue;
            for(Double d : gchisto.getChrHistogram(chr).retrieveRDBins()){
                sum += d;
                bins[d.intValue()] += 1;
            }
        }
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for(int i = 0; i < bins.length; i++){
            obs.add(i, bins[i]);
        }
        double testmean = sum / count;
        double teststdev = 0.0d;
        for(String chr : wins.getChrList()){
            for(Double d : gchisto.getChrHistogram(chr).retrieveRDBins()){
                teststdev += (double) Math.pow(d - testmean, 2.0d);
            }
        }
        teststdev /= (double) (count - 1);
        teststdev = Math.sqrt(teststdev);
        try{
            this.fitter = GaussianCurveFitter.create().withMaxIterations(50).withStartPoint(new double[]{maxvalue * 0.4 / teststdev,testmean, teststdev * 0.5});
        }catch(TooManyIterationsException ex){
            log.log(Level.WARNING, "Too many iterations! Using regular mean and stdev.");
            this.mean = testmean;
            this.stdev = teststdev;
            return;
        }
        double[] parameters = fitter.fit(obs.toList());
        Double mincut = parameters[1] - 2.0d * parameters[2];
        Double maxcut = parameters[1] + 2.0d * parameters[2];
        
        if(maxcut - mincut < 3 || mincut < 0 || maxcut < 0){
            log.log(Level.WARNING, "Gaussian fitting calculation had " + mincut + " and " + maxcut + "! Not fitting values");
            this.mean = parameters[1];
            this.stdev = parameters[2];
            return;
        }
        
        List<Double> tempvals = new ArrayList<>();
        for(int i = mincut.intValue(); i < maxcut.intValue(); i++){
            tempvals.add(bins[i]);
        }
        
        obs = new WeightedObservedPoints();
        for(int i = mincut.intValue(); i < maxcut.intValue(); i++){
        obs.add(i, bins[i]);
        }
        try{
            this.fitter = GaussianCurveFitter.create().withMaxIterations(50).withStartPoint(new double[]{maxvalue * 0.4 / teststdev,testmean, teststdev * 0.5});
        }catch(TooManyIterationsException ex){
            log.log(Level.WARNING, "Too many iterations! Using previously generated mean and stdev.");
            return;
        }
        double[] par = fitter.fit(obs.toList());
        this.mean = par[1];
        this.stdev = par[2];
    }
    
    public void CalculateMeanStdev(double[] values){
        Long maxvalue = 0l;
        if(values.length < 3){
            log.log(Level.WARNING, "Mean/Stdev calculation had " + values.length + "initial observations! Doing simple estimate");
            this.mean = StdevAvg.DoubleAvg(values);
            this.stdev = StdevAvg.stdevDBL(this.mean, values);
            return;
        }
        for(Double d : values){
            if(Math.round(d) > maxvalue)
                maxvalue = d.longValue();
        }
        
        if(maxvalue.intValue() <= 3){
            log.log(Level.WARNING, "Chromosome had very little signal! Max signal value: " + maxvalue.intValue());
            this.cantUse = true;
            return;
        }
        
        Double[] bins = new Double[maxvalue.intValue() + 1];
        java.util.Arrays.fill(bins, 0.0d);
        for(Double d : values){
            bins[d.intValue()] += 1;
        }
        double testmean = StdevAvg.DoubleAvg(values);
        double teststdev = StdevAvg.stdevDBL(testmean, values);
        
        log.log(Level.FINEST, "[GMSTDEV] testmean: " + testmean + " teststdev: " + teststdev);
        
        try{
            this.fitter = GaussianCurveFitter.create().withMaxIterations(50).withStartPoint(new double[]{maxvalue * 0.4 / teststdev,testmean, teststdev * 0.5});
        }catch(TooManyIterationsException ex){
            log.log(Level.WARNING, "Too many iterations for local fitter! Using regular mean and stdev.");
            this.mean = testmean;
            this.stdev = teststdev;
            return;
        }
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for(int i = 0; i < bins.length && i < testmean * 3; i++){
            obs.add(i, bins[i]);
        }
        
        double[] parameters = fitter.fit(obs.toList());
        
        this.mean = parameters[1];
        this.stdev = parameters[2];
        
        
        Double mincut = parameters[1] - 2.0d * parameters[2];
        Double maxcut = parameters[1] + 2.0d * parameters[2];
        
        log.log(Level.FINEST, "[GMSTDEV] initial mean: " + this.mean + " initial stdev: " + this.stdev + " mincut, maxcut: ( " + mincut + "," + maxcut + " )");
        
        if(maxcut - mincut < 3 || mincut < 0 || maxcut < 0){
            log.log(Level.WARNING, "Mean/Stdev calculation had " + mincut + " and " + maxcut + "! Not cropping values");
            return;
        }
        //java.util.Arrays.fill(bins, 0.0d);
        obs = new WeightedObservedPoints();
        for(int i = mincut.intValue(); i < maxcut.intValue(); i++){
            obs.add(i, bins[i]);
        }
        double[] par = fitter.fit(obs.toList());
        this.mean = par[1];
        this.stdev = par[2];
    }
    
    private void gaussianFit(double[] values, double constant, double mean, double sigma){
        double sum = 0.0d;
        double ss = 0.0d;
        for(double d : values){
            double v = (d - mean) / sigma;
            double x = constant * Math.exp(-0.5 * v * v);
            sum += x;
            
            ss += x * x;
        }
        
        this.fmean = sum / values.length;
        this.fstdev = Math.sqrt((ss / (double) values.length) - (this.fmean * this.fmean));
    }
    
    private void gaussianFitRange(double[] values, double constant){
        double min = this.fmean - 2.0d * this.fstdev;
        double max = this.fmean + 2.0d * this.fstdev;
        
        double sum = 0.0d;
        double ss = 0.0d;
        int n =0;
        for(double d : values){
            double v = (d - this.fmean) / this.fstdev;
            double x = constant * Math.exp(-0.5 * v * v);
            if(x >= min && x <= max){  
                sum += x;            
                ss += x * x;
                n++;
            }
        }
        
        this.mean = sum / n;
        this.stdev = Math.sqrt((ss / (double) n) - (this.mean * this.mean));
    }
    
    public double getMean(){
        return this.mean;
    }
    public double getStdev(){
        return this.stdev;
    }
}
