/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package FastaUtils;

import DataUtils.TempHistogram;
import java.nio.file.Path;

/**
 *
 * @author Derek.Bickhart
 */
public class GCHistogram extends TempHistogram<Double>{

    public GCHistogram(String chr, Path tmpdir) {
        super(chr, tmpdir);
    }

    @Override
    public void addHistogram(String chr, int start, int end) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void addHistogram(String chr, int start, int end, Double score) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void writeToTemp() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void readTemp() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
