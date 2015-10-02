/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils;

import file.BedAbstract;
import implement.BedSimple;

/**
 *
 * @author desktop
 */
public class BedStats extends BedSimple{
    private double copynum;
    private double eValue;
    private double gausEvalue;
    public BedStats(String c, int s, int e, String name, double copynum, double eValue, double gausEvalue) {
        super(c, s, e, name);
        this.copynum = copynum;
        this.eValue = eValue;
        this.gausEvalue = gausEvalue;
    }
    
    public String getOutputString(){
        StringBuilder str = new StringBuilder();
        str.append(chr).append("\t").append(start).append("\t").append(end);
        str.append("\t").append(name).append("\t").append(copynum).append("\t").append(eValue).append("\t");
        str.append(gausEvalue).append(System.lineSeparator());
        return str.toString();
    }
    
    @Override
    public int compareTo(BedAbstract o){
        if(!Chr().equals(o.Chr()))
            return Chr().hashCode() - o.Chr().hashCode();
        return Start() - o.Start();
    }
}
