/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author bickhart
 */
public class StdevAvg {
    public static double IntAvg(List<Integer> sum){
        if(sum.isEmpty()){
            return 0.0d;
        }
        double s = 0.0d;
        for(int d : sum){
            s += d;
        }
        return s / (double) sum.size();
    } 
    public static double DoubleAvg(double[] sum){
        if(sum.length == 0){
            return 0.0d;
        }
        double s = 0.0d;
        for(int x = 0; x < sum.length; x++){
            s += sum[x];
        }
        return s / sum.length;
    }
    public static double DoubleAvg(List<Double> sum){
        if(sum.isEmpty()){
            return 0.0d;
        }
        double s = 0.0d;
        for(double d : sum){
            s += d;
        }
        return s / sum.size();
    }
    public static float FloatAvg (float[] sum){
        if(sum.length == 0){
            return 0.0f;
        }
        float s = 0.0f;
        for(int x = 0; x < sum.length; x++){
            s += sum[x];
        }
        return s / sum.length;
    }
    
    public static double convertFltAvg(List<Float> sum){
        if(sum.size() == 0){
            return 0.0d;
        }
        double d = 0.0d;
        for(int x = 0; x < sum.size(); x++){
            d += (double) sum.get(x);
        }
        return d / (double) sum.size();
    }
    
    public static double stdevFlt(List<Float> sum){
        if(sum.size() == 0){
            return 0;
        }
        double mean = convertFltAvg(sum);
        double dev = 0.0d;
        for(int x = 0; x < sum.size(); x++){
            dev += (double) Math.pow(sum.get(x) - mean, 2.0d);
        }
        double variance = dev / (double) (sum.size() - 1);
        return Math.sqrt(variance);
    }
    
    public static double stdevFlt(double avg, List<Float> sum){
        if(sum.isEmpty() || sum.size() == 1){
            return 0;
        }        
        double dev = 0.0d;
        for(int x = 0; x < sum.size(); x++){
            dev += (double) Math.pow(sum.get(x) - avg, 2.0d);
        }
        double variance = dev / (double) (sum.size() - 1);
        return Math.sqrt(variance);
    }
    
    public static double stdevInt(double avg, List<Integer> sum){
        if(sum.isEmpty() || sum.size() == 1){
            return 0;
        }        
        double dev = 0.0d;
        for(int x = 0; x < sum.size(); x++){
            dev += (double) Math.pow(sum.get(x) - avg, 2.0d);
        }
        double variance = dev / (double) (sum.size() - 1);
        return Math.sqrt(variance);
    }
    
    public static double stdevDBL(double avg, List<Double> sum){
        if(sum.isEmpty() || sum.size() == 1){
            return 0;
        }        
        double dev = 0.0d;
        for(double d : sum){
            dev += (double) Math.pow(d - avg, 2.0d);
        }
        double variance = dev / (double) (sum.size() - 1);
        return Math.sqrt(variance);
    }
    
    public static double stdevDBL(double avg, double[] sum){
        if(sum.length <= 1){
            return 0;
        }        
        double dev = 0.0d;
        for(double d : sum){
            dev += (double) Math.pow(d - avg, 2.0d);
        }
        double variance = dev / (double) (sum.length - 1);
        return Math.sqrt(variance);
    }
    
    public static double getRangeAverage(double[] values, int start, int end){
        // End index inclusive
        int len = end - start;
        if(len <= 0 || values.length == 0)
            return 0.0d;
        
        double sum = 0.0d;
        for(int i = start; i <= end; i++){
            sum += values[i];
        }
        return sum / len;
    }
    
    public static double getRangeVariance(double[] values, double average, int start, int end){
        // End index inclusive
        int len = end - start;
        if(len <= 0 || values.length == 0)
            return 0.0d;
        
        double ss = 0.0d;
        for(int i = start; i <= end; i++){
            ss += Math.pow(values[i] - average, 2.0d);
        }        
        return ss / (end - start);
    }
}
