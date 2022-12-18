package riemann;

import static org.junit.Assert.*;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import math.GSeries;

public class RiemannTest {
    static final int initialPadding = 40;

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
    
    @Test
    public void testCrossProduct() throws IOException {
        double oldZeta = Double.NEGATIVE_INFINITY;
        double crossSum = 0;
        int sample = 0;
        double[] gramSum = {0,0};
        double t0 = StaticMethods.gramE12[0][0];
        for (int i = 0; i < StaticMethods.gramE12.length; i++) {
            double[] limits = StaticMethods.gramE12[i];
            GSeries gSeries = StaticMethods.getSavedGSeries(limits[0], BigDecimal.valueOf(1.0E12));
            double tIncr = gSeries.spacing;
            for (int j = 0; j < 30100; j++) {
                double zeta = gSeries.evaluateZeta(t0, initialPadding);
                gramSum[sample%2] += zeta;
                if(oldZeta != Double.NEGATIVE_INFINITY) {
                    crossSum += oldZeta*zeta;
                    sample++;
                }
                oldZeta = zeta;
                t0 += tIncr;
                if(t0 > limits[1]) {
                    break;
                }
            }
            double[] gramMean = {gramSum[0]/sample, gramSum[1]/sample};
            System.out.println("mean crossSum " + crossSum/sample
                + " sample " + sample + " " + Arrays.toString(gramMean));
        }
        System.out.println("mean crossSum " + crossSum/sample
         + " sample " + sample);
    }
    
    @Test
    public void testPoly7() throws IOException {
        /*
        deviation 1.0025324856721696E-4 Bad Poly7{
        a=246.09380502666446, b=246.32312342197378, c=246.7424183211807,
        d0=-17.718867338654096, d1=27.697678463431405, d2=-10.289947584846136,
        t0=-800.9150900144672, t1=2995.8938595958193, t2=-139.12408210913475,
         m0=-1.594159242026056, m1=4.010423055063137, offset=246.32312342197378}
246.09380502666446,
246.32312342197378,
246.7424183211807,

-17.718867338654096,
27.697678463431405,
-10.289947584846136,

-1.594159242026056,
4.010423055063137,
246.32312342197378}

deviation NaN Bad Poly7{
366.6474495506842,
366.9272070782925,
367.1694973548008,
9.49029728124495,
-14.602482306322326,
6.709988520987169,

1.3190132125782246,
-0.9281998769066759,
366.9272070782925
}

         */
        Poly7.epsilon = 1.0E-6;
        Poly7.derepsilon = 1.0E-6;
        double[] limits = StaticMethods.gramE12[0];
        double t = limits[0];
        GSeries gSeries = StaticMethods.getSavedGSeries(limits[0], BigDecimal.valueOf(1.0E12));
        System.out.println("zeta " +
            gSeries.evaluateZeta(t, initialPadding) +
            " t  " + t );
        double[][] zeroInfo = {
            {244.15890691298068396, -20.007604626096071598,  -1.232146174810101691},
            {244.367502584863394599,  19.343950349024609636,  1.554959200487025184},
            {244.588579452072075626,  -27.175877175067416402,  -3.420984658449017779},
        };
        Poly7 poly7 = new Poly7(
            zeroInfo[0][0], zeroInfo[1][0], zeroInfo[2][0],
            zeroInfo[0][1], zeroInfo[1][1], zeroInfo[2][1]
            
        );
        poly7.setExtrema(zeroInfo[0][2], zeroInfo[1][2], zeroInfo[1][0]);
        poly7.setTermValues();
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        System.out.println("poly7term " + poly7.poly7term.A
            + " " + poly7.poly7term.B);
        System.out.println("============= " );
        
        double incr  = gSeries.spacing/2;
        t += 3*incr;
        for (int i = 0; i < 3; i++) {
            t += incr;
            System.out.println("zeta " +
                gSeries.evaluateZeta(t, initialPadding) +
                " t  " + t );
            System.out.println("poly " +
                poly7.eval(t)
                );
        }
    }
    
    @Test
    public void badCase1() {
        Poly7 poly7 = new Poly7(
            366.6474495506842,
            366.9272070782925,
            367.1694973548008,
            9.49029728124495,
            -14.602482306322326,
            6.709988520987169
        );
        poly7.tabulate(
            366.927207078292,
            367.1694973548008,
            20
        );
        poly7.setExtrema(
            1.3190132125782246,
            -0.9281998769066759,
            366.9272070782925
        );
        double deviation = poly7.setTermValues();
        System.out.println("deviation " + deviation);
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        System.out.println("============= " );
        
    }
    
    @Test
    public void badCase2() {
        Poly7.epsilon = 1.0E-6;
        Poly7.derepsilon = 1.0E-6;
        /*
[243.8749480149, 244.15890691298068, 244.3675025848634]
[17.619276585379914, -20.00760462609607, 19.34395034902461]
Rosser.extrema [1.9266754104451154, -1.2321461748101017, 1.5549592004870252]
        
         */
        Poly7 poly7 = new Poly7(
            244.3675025848634, 244.58857945207208, 244.92059950582586,
            19.34395034902461, -27.175877175067416, 23.85164367971759
        );
        poly7.setExtrema(
            1.5549592004870252,3.420984658449018, 244.58857945207208
        );
        double deviation = poly7.setTermValues();
        System.out.println("deviation " + deviation);
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        System.out.println("============= " );
        
    }
    
    @Test
    public void testSavedE12() throws IOException {
        for (int i = 0; i < StaticMethods.gramE12.length; i++) {
            double[] limits = StaticMethods.gramE12[i];
            System.out.println(Arrays.toString(limits));
            System.out.println("idx  " +
                (limits[0]-244.02115917156451839965694310614387)/0.24359904590398668 +
                " idx at zero " + (limits[1]-244.02115917156451839965694310614387)/0.24359904590398668);
            GSeries gSeries = StaticMethods.getSavedGSeries(limits[0], BigDecimal.valueOf(1.0E12));
            System.out.println("zeta at i  " + i +   " " +
                gSeries.evaluateZeta(limits[0], initialPadding) +
                " t  " +
                limits[0] + ", " + (limits[0]  + 7307.97133775572));
            if(i < StaticMethods.gramE12.length-1) {
                System.out.println("zeta at next after " + i +  " "
                    + gSeries.evaluateZeta(StaticMethods.gramE12[i+1][0], initialPadding));
            }
        }
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
        System.out.println("incr " + incr);
        
        BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
        double predictedSqrtArg1 = sqrtArg1.doubleValue();
        double delSqrtArg1 = incr/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
        double cos = -1;
        // math.MoreGSeriesTest.test1E12
        System.out.println("n-3945951431270L,  zeta");
        double txx = begin;
        for (int i = 1; i <= R; i++) {
            predictedSqrtArg1 += delSqrtArg1;
            double rotatedSum = 2*( cos*fAtBeta[i-1][0]);
            double correction = GSeries.correction( predictedSqrtArg1);
            double zeta = rotatedSum + correction;
            System.out.println(i + ", " + zeta+ ", " + (txx));
            txx += incr;
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
