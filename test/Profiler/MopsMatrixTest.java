/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Profiler;

import java.util.List;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author dbickhart
 */
public class MopsMatrixTest {
    
    public MopsMatrixTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }

    /**
     * Test of MopsAlgorithm method, of class MopsMatrix.
     */
    @Test
    public void testMopsAlgorithm() {
        System.out.println("MopsAlgorithm");
        int normStateCN = 0;
        int cycNum = 0;
        MopsMatrix instance = null;
        List<INISeg> expResult = null;
        List<INISeg> result = instance.MopsAlgorithm(normStateCN, cycNum);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
    
}
