/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jarms.modes;

import GetCmdOpt.ArrayModeCmdLineParser;
import java.util.logging.Logger;

/**
 *
 * @author desktop
 */
public class ProfileMode {
    private static final Logger log = Logger.getLogger(ProfileMode.class.getName());
    private final ArrayModeCmdLineParser cmd;
    
    public ProfileMode(ArrayModeCmdLineParser cmd){
        this.cmd = cmd;
    }
}
