package riemann;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import math.GSeries;

public class RiemannTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
    /**
     * Test method for Riemann zeta at gram points 1.0E12.
     */
    @Test
    public void test1E12() {
        int k0 = 1, k1=398942;
        Gram.initLogVals(k1);
        int R = 10;
        long init= System.currentTimeMillis();
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        //gram [GramInfo [grampt=1000000000243.77756012466052947405878015472510, idx=3945951431271], 
        //      GramInfo [grampt=1000000000244.02115917156451839965694310614387, idx=3945951431272]]
        //0.24359904690398881 tincr 0.24359904690399015
        double begin = 243.77756012466052947405878015472510;
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
        BigDecimal tval = new BigDecimal(begin, Gram.mc).add(
              offset, Gram.mc);
        double[][] fAtBeta = GSeries.fSeries(k0, k1, incr, R, tval);

        long end = System.currentTimeMillis();
        System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
        
        BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
        double predictedSqrtArg1 = sqrtArg1.doubleValue();
        double delSqrtArg1 = incr/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
        double cos = -1;
        System.out.println("n-3945951431270L,  zeta");
        for (int i = 1; i <= R; i++) {
            predictedSqrtArg1 += delSqrtArg1;
            double rotatedSum = 2*( cos*fAtBeta[i-1][0]);
            double correction = GSeries.correction( predictedSqrtArg1);
            double zeta = rotatedSum + correction;
            System.out.println(i + ", " + zeta);
            cos = -cos;
        }
        
    }
    
	@Test
	public void testRiemannDoubleLong() {
		double[] args = {1.87383225, 2.13054990,        
		        2.38726754,   2.64398519, 1.8475231278        };
		double[] vals = {4.023671,42.85232072, -4.52434707, 
				-2.84388877, 0};
		long init = 0;
		for (int i = 0; i < args.length; i++) {
			double zeta = Riemann.riemann(args[i], 267653395647L);
			System.out.println(args[i] + ":" +zeta + ":" + Math.abs(zeta-vals[i]));
			assertTrue(i + ":" + Math.abs(zeta-vals[i]), Math.abs(zeta-vals[i])<0.0000005);
			if(i==0){init = System.currentTimeMillis();}
		}
		long end = System.currentTimeMillis();
		System.out.println("N " + Riemann.oldN + " calc for " + (args.length-1) + " " + (end - init) + "ms");
	}

}
