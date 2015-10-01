/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

import DataUtils.BinCoords;

/**
 *
 * @author desktop
 */
public class EvalueTools {
    
    public static boolean AdjustToEvalue(double mean, double sd, int genomeSize, TDistributionFunction tdist, double[] rd, BinCoords b, double eval){
        final int MAX_STEPS = 1000;
        int s[] = new int[MAX_STEPS], e[] = new int[MAX_STEPS];
        int step_count = 0;
        s[step_count] = b.start;
        e[step_count] = b.end;
        step_count++;
        int ll_ind = 1,lr_ind = 2,rl_ind = 3,rr_ind = 4;

        double val = GetEValue(mean,sd, genomeSize, tdist, rd,b.start,b.end);
        while (b.end > b.start + 1 && val > eval && step_count < MAX_STEPS) {
          double best_e = 1e+10, tmp = 0;
          int best_index = 0;
          if (b.start - 1 >= 0 && // Left left
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,b.start - 1,b.end)) < best_e) {
            best_e = tmp;
            best_index = ll_ind;
          }
          if (b.start + 1 < rd.length && // Left right
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,b.start + 1,b.end)) < best_e) {
            best_e = tmp;
            best_index = lr_ind;
          }
          if (b.end - 1 >= 0 && // Right left
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,b.start,b.end - 1)) < best_e) {
            best_e = tmp;
            best_index = rl_ind;
          }
          if (b.end + 1 < rd.length && // Right right
              (tmp = GetEValue(mean,sd,genomeSize, tdist,rd,b.start,b.end + 1)) < best_e) {
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

        if (b.end > b.start && val <= eval) return true;
        return false;
    }
    
    public static void GetGaussianEValue(){
        
    }
    
    public static double GetEValue(double mean, double sd, int genomeSize, TDistributionFunction tdist, double[] rd, int start, int end){
        int n = end - start + 1;
        double aver = 0,s = 0, over_n = 1./n;
        for (int b = start;b <= end;b++) {
          aver  += rd[b];
          s += rd[b]*rd[b];
        }
        aver *= over_n;
        s = Math.sqrt(s*over_n - aver*aver);

        if (s == 0) s = sd * Math.sqrt(aver/mean);
        if (s == 0) s = 1;
        
        double x = (aver - mean)* Math.sqrt(n)/s;
        double ret = tdist.getDensityValue(n, x);

        if (x > 0) ret = 1 - ret;
        
        //CHECK: originally, this was "*" the inverse of the bin count and "*" the length of the bin stretch
        ret *= genomeSize / (rd.length * (end - start + 1));

        return ret;
    }
    
}
