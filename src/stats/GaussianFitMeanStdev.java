/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import DataUtils.WindowPlan;
import HistogramUtils.ChrHistogramFactory;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

/**
 *
 * @author Derek.Bickhart
 */
public class GaussianFitMeanStdev {
    GaussianCurveFitter fitter = GaussianCurveFitter.create();
    private double fmean = 0.0d;
    private double fstdev = 0.0d;
    private double mean;
    private double stdev;
    
    public void CalculateGlobalMeanStdev(ChrHistogramFactory gchisto, WindowPlan wins){
        Long maxvalue = 0l;
        for(String chr : wins.getChrList()){
            if(!gchisto.hasChrHistogram(chr))
                continue;
            for(Double d : gchisto.getChrHistogram(chr).retrieveRDBins())
                if(d.longValue() > maxvalue)
                    maxvalue = d.longValue();
        }
        
        Double[] bins = new Double[maxvalue.intValue() + 1];
        java.util.Arrays.fill(bins, 0.0d);
        for(String chr : wins.getChrList()){
            if(!gchisto.hasChrHistogram(chr))
                continue;
            for(Double d : gchisto.getChrHistogram(chr).retrieveRDBins())
                bins[d.intValue()] += 1;
        }
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for(int i = 0; i < bins.length; i++){
            obs.add(i, bins[i]);
        }
        double[] parameters = fitter.fit(obs.toList());
        Double mincut = parameters[1] - 2.0d * parameters[2];
        Double maxcut = parameters[1] + 2.0d * parameters[2];
        obs = new WeightedObservedPoints();
        for(int i = mincut.intValue(); i < maxcut.intValue(); i++){
            obs.add(i, bins[i]);
        }
        
        double[] par = fitter.fit(obs.toList());
        this.mean = par[1];
        this.stdev = par[2];
    }
    
    public void CalculateMeanStdev(double[] values){
        Long maxvalue = 0l;
        for(Double d : values){
            if(Math.round(d) > maxvalue)
                maxvalue = d.longValue();
        }
        
        Double[] bins = new Double[maxvalue.intValue() + 1];
        java.util.Arrays.fill(bins, 0.0d);
        for(Double d : values){
            bins[d.intValue()] += 1;
        }
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for(int i = 0; i < bins.length; i++){
            obs.add(i, bins[i]);
        }
        double[] parameters = fitter.fit(obs.toList());
        
        this.mean = parameters[1];
        this.stdev = parameters[2];
        
        
        Double mincut = parameters[1] - 2.0d * parameters[2];
        Double maxcut = parameters[1] + 2.0d * parameters[2];
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
