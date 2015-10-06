/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package DataUtils;

/**
 *
 * @author desktop
 */
public class BinCoords {
    public int start = 0;
    public int end = 0;
    public boolean useable = false;
    
    public boolean isNormalInterval(){
        return this.end > this.start;
    }
    
    public int getLength(){
        return this.end - this.start;
    }
}
