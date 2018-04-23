/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Profiler;

import DataUtils.WindowPlan;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import org.apache.commons.math3.special.Gamma;
import stats.StdevAvg;

/**
 *
 * @author dbickhart
 */
public class MopsMatrix {
    private final double[][] chrRDMat;
    private final double[][][] alpha_ik;
    private final List<String> sampleOrder;
    private final String chr;
    private final double[] I;
    private final Map<String, Double> cov = new HashMap<>();
    
    public MopsMatrix(SampleList data, WindowPlan wins, List<String> sampleOrder, String chr, final double[] I) throws Exception{
        this.sampleOrder = sampleOrder;
        this.chr = chr;
        this.I = I;
        
        this.chrRDMat = new double[wins.getNumBins(chr)][sampleOrder.size() -1];
        this.alpha_ik = new double[wins.getNumBins(chr)][(sampleOrder.size() -1)][I.length - 1]; 
        int arrayIdx = 0;
        for(String s : this.sampleOrder){
            double[] rdVals = data.getRDValues(s, chr);
            // Calculate values for coverage estimates
            final double average = StdevAvg.DoubleAvg(rdVals);
            final double stDev = StdevAvg.stdevDBL(average, rdVals);
            // Simple exclusion of outliers based on Stdevs
            double[] filtered = Arrays.stream(rdVals)
                            .filter(d -> d > average - (3 * stDev) && d < average + (3 * stDev))
                            .toArray();
            final double updatedAvg = StdevAvg.DoubleAvg(filtered);
            double median = StdevAvg.getMedian(filtered);
            
            if(median == 0)
                throw new Exception("Sample: " + s + " median coverage was found to be Zero for chr: " + chr);
            
            this.cov.put(s, updatedAvg / median);
            
            // Filling the RD array with uncorrected RD values
            for(int x = 0; x < wins.getNumBins(chr); x++)
                this.chrRDMat[x][arrayIdx] = rdVals[x];
            
            arrayIdx++;
        }
    }
    
    public Map<String, List<Double[]>> MopsAlgorithm(int normStateCN, int cycNum){
        // The first double is the I/NI and the second is the signed I/NI
        Map<String, List<Double[]>> INIvals = new ConcurrentHashMap<>(sampleOrder.size());
        
        // Setting the initial alpha values
        double[] alphaInit = new double[this.I.length - 1];
        for(int x = 0; x < this.I.length; x++)
            alphaInit[x] = (x == normStateCN)? 0.6 : 0.05;
        
        double ainitSum = Arrays.stream(alphaInit).sum();
        for(int x = 0; x < alphaInit.length; x++)
            alphaInit[x] /= ainitSum;
        
        // Setting the alpha prior value
        double[] alphaPrior = new double[this.I.length -1];
        for(int x = 0; x < alphaPrior.length; x++)
            alphaPrior[x] = (x == normStateCN)? 1 : 0;
        
        // Setting the lambda initial state
        double[] covVec = new double[this.sampleOrder.size() - 1];
        int i = 0;
        for(String s : this.sampleOrder)
            covVec[i++] = this.cov.get(s);
        double[][] lambdaEst = divideMatVec(this.chrRDMat, covVec);
        double[][] lambdaInit = divideMatVec(lambdaEst, this.I);
        
        double[][] lgamma = new double[this.chrRDMat.length][this.chrRDMat[0].length];
        for(int x = 0; x < this.chrRDMat.length; x++)
            for(int y = 0; y < this.chrRDMat[0].length; y++){
                lgamma[x][y] =  Gamma.logGamma(this.chrRDMat[x][y]); // NOTE: this function does not return gamma values for zero!
            }
        
        
        
        return INIvals;
    }
    
    private double[][] multiplyMatVec(double[][] M, double[] V){
        if(M.length != V.length)
            throw new IllegalArgumentException("M:Rows did not match V:Rows!");
        
        double[][] C = new double[M.length][M[0].length];
        /*for(int i = 0; i < M.length; i++)
        for(int j = 0; j < M[0].length; j++)
        C[i][j] = 0.0;*/
        
        for(int i = 0; i < M.length; i++)
            for(int j = 0; j < M[0].length; j++)
                C[i][j] = M[i][j] * V[i];
        
        return C;
    } 
    
    private double[][] divideMatVec(double[][] M, double[] V){
        if(M.length != V.length)
            throw new IllegalArgumentException("M:Rows did not match V:Rows!");
        
        double[][] C = new double[M.length][M[0].length];
        /*for(int i = 0; i < M.length; i++)
        for(int j = 0; j < M[0].length; j++)
        C[i][j] = 0.0;*/
        
        for(int i = 0; i < M.length; i++)
            for(int j = 0; j < M[0].length; j++)
                C[i][j] = (V[i] > 0)? M[i][j] / V[i] : 0;
        
        return C;
    } 
    private Double[][] multiplicar(Double[][] A, Double[][] B) {

        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;

        if (aColumns != bRows) {
            throw new IllegalArgumentException("A:Rows: " + aColumns + " did not match B:Columns " + bRows + ".");
        }

        Double[][] C = new Double[aRows][bColumns];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                C[i][j] = 0.00000;
            }
        }

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }
}
