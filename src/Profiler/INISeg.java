/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Profiler;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * @author dbickhart
 */
public class INISeg {

    public final double INI;
    private final Map<String, Double> SINI = new ConcurrentHashMap<>();

    public INISeg(double INI, double[] SINI, List<String> samples){
        this.INI = INI;
        for(int x = 0; x < samples.size(); x++)
            this.SINI.put(samples.get(x), SINI[x]);
    }
    public INISeg(double INI, List<Double> SINI, List<String> samples){
        this.INI = INI;
        for(int x = 0; x < samples.size(); x++)
            this.SINI.put(samples.get(x), SINI.get(x));
    }

    public Double getSINI(String sample){
        return this.SINI.get(sample);
    }

}
