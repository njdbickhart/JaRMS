/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;


import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;

/**
 * The T test to determine if a partition is significantly different from surrounding partitions
 * @author desktop
 */
public class TDistributionFunction {
    private static final Logger log = Logger.getLogger(TDistributionFunction.class.getName());
    
    private final Map<Integer, TDistribution> functions = new ConcurrentHashMap<>();
    private final int winsize;
    
    
    public TDistributionFunction(int winsize){
        // This is the constant window size for the bins
        this.winsize = winsize;
    }
    
    public synchronized double getDensityValue(int n, double value){
        if(n <= 1)
            return 1.0d;
        try{
            if(!this.functions.containsKey(n))
                functions.put(n, new TDistribution(n));
        }catch(NotStrictlyPositiveException ex){
            log.log(Level.SEVERE, "Error attempting to create TDist function for: " + n, ex);
        }
        return this.functions.get(n).cumulativeProbability(value);
    }
    
    public double TestOneRegion(double value, double mean, double sigma, int n){
        double x = (value - mean) / Math.sqrt(sigma) * Math.sqrt(n);
        double p = this.getDensityValue(n -1, x); 
        if (x > 0) p = 1 - p;
        return p;
    }
    
    public double TestTwoRegions(double m1, double s1, int n1, double m2, double s2, int n2, double scale){
        if (s1 == 0) s1 = 1;
        if (s2 == 0) s2 = 1;
        double tmp1 = s1*(1.0d / n1), tmp2 = s2*(1.0d / n2);
        double s = Math.sqrt(tmp1 + tmp2);
        double t = (m1 - m2)/s;
        double tmp = (tmp1 + tmp2)*(tmp1 + tmp2)*(n1 - 1)*(n2 - 1);
        tmp /= tmp1*tmp1*(n2 - 1) + tmp2*tmp2*(n1 - 1);
        int ndf = (int) Math.round(tmp + 0.5);
        
        if(ndf <= 1){
            System.out.println("Error finding proper T function!");
        }
        
        double ret = this.getDensityValue(ndf, t);
        if (t > 0 && ret != 1d) 
            ret = 1 - ret;

        ret *= scale * (1.0d / this.winsize) * (1.0d /(n1 + n2));

        return ret;
    }
}
