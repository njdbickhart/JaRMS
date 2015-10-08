/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import DataUtils.BinCoords;
import org.apache.commons.math3.special.Erf;

/**
 *
 * @author desktop
 */
public class EvalueTools {
    
    public static BinCoords AdjustToEvalue(double mean, double sd, long genomeSize, TDistributionFunction tdist, double[] rd, BinCoords b, int winsize, double eval){
        final int MAX_STEPS = 1000;
        b.useable = false;
        int s[] = new int[MAX_STEPS], e[] = new int[MAX_STEPS];
        int step_count = 0;
        s[step_count] = b.start;
        e[step_count] = b.end;
        step_count++;
        int ll_ind = 1,lr_ind = 2,rl_ind = 3,rr_ind = 4;

        double val = GetEValue(mean,sd, genomeSize, tdist, rd,winsize,b.start,b.end);
        while (b.end > b.start + 1 && val > eval && step_count < MAX_STEPS) {
          double best_e = 1e+10, tmp = 0;
          int best_index = 0;
          if (b.start - 1 >= 0 && // Left left
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,winsize,b.start - 1,b.end)) < best_e) {
            best_e = tmp;
            best_index = ll_ind;
          }
          if (b.start + 1 < rd.length && // Left right
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,winsize,b.start + 1,b.end)) < best_e) {
            best_e = tmp;
            best_index = lr_ind;
          }
          if (b.end - 1 >= 0 && // Right left
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,winsize,b.start,b.end - 1)) < best_e) {
            best_e = tmp;
            best_index = rl_ind;
          }
          if (b.end + 1 < rd.length && // Right right
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,winsize,b.start,b.end + 1)) < best_e) {
            best_e = tmp;
            best_index = rr_ind;
          }
          if (best_e > val) break; // Can't improve e-value
          if (best_index == ll_ind)      b.start -= 1;
          else if (best_index == lr_ind) b.start += 1;
          else if (best_index == rl_ind) b.end   -= 1;
          else if (best_index == rr_ind) b.end   += 1;
          val = best_e;

          for (int i = 0;i < step_count;i++)
            if (b.start == s[i] && b.end == e[i]) break; // Get into loop
          s[step_count] = b.start;
          e[step_count] = b.end;
          step_count++;
        }

        if (b.end > b.start && val <= eval) 
            b.useable = true;
        return b;
    }
    
    public static double GetGaussianEValue(double mean,double sigma,double rd[],
				int start,int end, long GenomeSize){
        // Calculate by deviation from gaussian
        double max = 0,min = 1e+10,av = 0;
        int n = end - start + 1;
        for (int i = start;i <= end;i++) {
          av += rd[i];
          if (rd[i] > max) max = rd[i];
          if (rd[i] < min) min = rd[i];
        }

        av /= n;

        double p;
        if (av < mean) {
          double x = (max - mean)/sigma*0.707;
          p = 0.5*(1 + Erf.erf(x));
        } else {
          double x = (min - mean)/sigma*0.707;
          p = 0.5*(1 - Erf.erf(x));
        }
        return GenomeSize * Math.pow(p, n);
    }
    
    public static double GetEValue(double mean, double sd, long genomeSize, TDistributionFunction tdist, double[] rd, int winsize, int start, int end){
        int n = end - start + 1;
        double aver = 0,s = 0, ss = 0, over_n = 1./n;
        for (int b = start;b <= end;b++) {
          aver  += rd[b];
          s += rd[b]*rd[b];
        }
        aver *= over_n;
        //NOTE: this is a test for estimating the stdev using an alternate method
        for(int b = start; b <= end; b++){
            ss += Math.pow(rd[b] - aver, 2.0d);
        }
        
        s = Math.sqrt(s*over_n - aver*aver);
        
        // NOTE: still part of my calculation test
        ss /= n;
        ss = Math.sqrt(ss);
        
        if (s == 0) s = sd * Math.sqrt(aver/mean);
        if (s == 0) s = 1;
        
        double x = (aver - mean)* Math.sqrt(n)/s;
        double ret = tdist.getDensityValue(n - 1, x);
        
        // NOTE: still part of my calculation test
        double testret = tdist.getDensityValue(n, ss);
        
        if (x > 0) ret = 1 - ret;
        
        
        ret *= (genomeSize * (1 - 0.10)) / (winsize * (end - start + 1));
        
        
        return ret;
    }
    
}
