/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package stats;

/**
 *
 * @author desktop
 */
public class EvalueTools {
    
    public static boolean AdjustToEvalue(double mean, double sd, double[] rd, BinCoords b, double eval){
        final int MAX_STEPS = 1000;
        int s0 = b.start, e0 = b.end;
        int s[MAX_STEPS],e[MAX_STEPS],step_count = 0;
        s[step_count] = b.start;
        e[step_count] = b.end;
        step_count++;
        int ll_ind = 1,lr_ind = 2,rl_ind = 3,rr_ind = 4;

        double val = getEValue(mean,sd,rd,b.start,b.end);
        while (b.end > b.start + 1 && val > eval && step_count < MAX_STEPS) {
          double best_e = 1e+10, tmp = 0;
          int best_index = 0;
          if (b.start - 1 >= 0 && // Left left
              (tmp = getEValue(mean,sd,rd,b.start - 1,b.end)) < best_e) {
            best_e = tmp;
            best_index = ll_ind;
          }
          if (b.start + 1 < n_bins && // Left right
              (tmp = getEValue(mean,sd,rd,b.start + 1,b.end)) < best_e) {
            best_e = tmp;
            best_index = lr_ind;
          }
          if (b.end - 1 >= 0 && // Right left
              (tmp = getEValue(mean,sd,rd,b.start,b.end - 1)) < best_e) {
            best_e = tmp;
            best_index = rl_ind;
          }
          if (b.end + 1 < n_bins && // Right right
              (tmp = getEValue(mean,sd,rd,b.start,b.end + 1)) < best_e) {
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
    
    public static double void GetEValue(){
        int n = end - start + 1;
        double aver = 0,s = 0, over_n = 1./n;
        for (int b = start;b <= end;b++) {
          aver  += rd[b];
          s += rd[b]*rd[b];
        }
        aver *= over_n;
        s = TMath::Sqrt(s*over_n - aver*aver);

        if (s == 0) s = sigma*TMath::Sqrt(aver/mean);
        if (s == 0) s = 1;

        TF1* cum = getTFunction(n - 1);
        double x = (aver - mean)*getSqrt(n)/s;
        double ret = cum->Eval(x);

        if (x > 0) ret = 1 - ret;
        ret *= GENOME_SIZE_NORMAL*inv_bin_size*getInverse(end - start + 1);

        return ret;
    }
    
    public class BinCoords{
        public int start = 0;
        public int end = 0;
    }
}
