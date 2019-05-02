package riemann;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import math.GSeries;
import riemann.Interpolate.GramOrMid;
import riemann.Interpolate.Poly3;
import riemann.Interpolate.Poly4;

public class InterpolateTest {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }
    
    public double eval(double x, Poly3 poly, double C){
        double mult = (x-poly.z0)*(x-poly.z1);
        mult = mult*mult/poly.denom;
        double ret = poly.eval1(x)+C*mult;
        return ret;
    }        

    @Test @Ignore
    public void testSlowGrowth() {
        final double z0 = 247.1270084390591, z1 = 247.25934052902014;
        final double d0 = -5.46251113952012, d1 = 4.963052066753181;
        final double max = 0.19014781349558818;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        double gram = 247.18794676831632;
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        System.out.println(nf.format(gram) 
                + ", " + nf.format(poly4.eval(gram))
                + ", der " + nf.format(poly4.der(gram))
                );
        assertEquals(-0.189428825055622, poly4.eval(gram), 0.00005);
        System.out.println();
    }
    
    @Test 
    public  void debug() {
        final double z0 = 1, z1 = 2;
        final double d0 = -1, d1 = 1;
        final double max = -0.25;
        Poly4 poly = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin = poly.der(poly.positionMax);
        System.out.println(poly.C + ", " + poly.positionMax + ", " 
        + poly.eval(poly.positionMax)+ ", der " + derAtMin);
        assertEquals(2.0d, poly.secondDerRHS(), 0.0000001);
        assertEquals(2.0d, poly.secondDer(1.75), 0.0000001);
        assertEquals(0d, poly.C, 0.0000001);
        assertEquals(0d, derAtMin, 0.0000001);
        int N = 11;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    );
        }
    }

    @Test 
    public void testE28() {
        final double z0 = 66093.28515358719, z1 = 66093.56881530148;
        final double d0 = 1866.9381635238121, d1 = -1323.4908633876328;
        final double max = 392.08238609990934;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        System.out.println(" positionMax " + poly4.positionMax
               +  ", der " + poly4.der(poly4.positionMax)
               + ", eval " + poly4.eval(poly4.positionMax));
        assertEquals(max, poly4.eval(poly4.positionMax),0.000001);
        assertEquals(0d, poly4.der(poly4.positionMax),0.000001);
    }
    

    @Test
    public void testInterpolate() {
        // x^4-16 - (x-2)*(625-16)/3 
        // 4x^3-303 :4.231173938822708
        // 12x 
        double z0 = 2, z1 = 5;
        double d0 = -271, d1 = 197;
        double xmin = 4.231173938822708;
        double max = (xmin*xmin-4)*(xmin*xmin+4)- (xmin-2)*303;
        Poly4 poly = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin = poly.der(poly.positionMax);
        System.out.println(poly.C + ", " + poly.positionMax + ", " 
        + poly.eval(poly.positionMax)+ ", der " + derAtMin);
        assertEquals(max, poly.eval(poly.positionMax),0.000001);
        assertEquals(0d, Math.abs(derAtMin), 0.0000001);
    }
    
    @Test @Ignore("Use E12")
    public void testReadItems(   ) throws Exception {
    	//run math.GSeriesTest.test1E12F() first
        double begin= Interpolate.baseLimit + 
        		(Interpolate.noffset-Interpolate.correction)* (Interpolate.gramIncr);
        GSeries gSeries1 = new GSeries(1, 0, Interpolate.offset, begin, Interpolate.gramIncr);
        
        Interpolate.zetaCorrection1 = GSeries.correction( gSeries1.basesqrtArg1);
        
        BigDecimal tvalsi = Interpolate.offset.add(BigDecimal.valueOf(begin), Gram.mc);
        BigDecimal gramIndex1 = Gram.theta(tvalsi, Gram.mc).divide(Gram.pi, Gram.mc);
        String[] line = Rosser.getParam("header").split("[-L]+");
        
        gramIndex1 = 
                gramIndex1.subtract(new BigDecimal(line[1]), Gram.mc);
        System.out.println( gSeries1.begin + ", zetaCorrection " + Interpolate.zetaCorrection1
                + ", gram index " + gramIndex1 + ", " + line[1]);
        
        int N = 20;
        Interpolate.zeroInput = Rosser.readZeros(
        		Interpolate.baseLimit, 
        		Interpolate.out, Interpolate.zeroIn, null);
        System.out.println(Arrays.toString(Interpolate.zeroInput.lastZero)  +
                ", " + Interpolate.baseLimit + ", " + 
        		Arrays.toString(Interpolate.zeroInput.nextValues));
        System.arraycopy(Interpolate.zeroInput.nextValues, 0, 
        		Interpolate.lastZeroSeen1, 0, Interpolate.zeroIn.length);
        //f at Mid. 0 -> real. 1 -> im
        Interpolate.imFmid = new double[N][2];
        //f at Gram. 0 -> real. 1 -> im
        Interpolate.fAtBeta = new double[N][2];
        Interpolate.gramDer = new double[N][2];
        double[] zetaMidMean = {0, 0};
        double[] zetaGramMean = {0, 0};
        double[] zetaMid = new double[N];
        double[] zetaGram = new double[N];
        int idx = 0;
        while (idx < N  ) {
        	//this n is actually n-1!!!
        	//idx = 0, n = 3 in zetaE12.csv
            int nprime = idx + Interpolate.noffset;
            double upperLimit = Interpolate.baseLimit + 
            		(nprime-Interpolate.correction-1)* (Interpolate.gramIncr);
            // populate fAtBeta,  zetaGramMean
            //factor of 2 in correction
            zetaGram[idx] =  2*Interpolate.getZetaEstimate(nprime, idx, upperLimit, 
            		zetaGramMean,	Interpolate.fAtBeta, GramOrMid.GRAM);

            
            upperLimit += Interpolate.gramIncr/2;
            zetaMid[idx] = 2*Interpolate.getZetaEstimate(nprime, idx, upperLimit, zetaMidMean,
            		Interpolate.imFmid, GramOrMid.MID);
            idx++;
        }
        System.out.println( "breaks: " + Interpolate.breaks);
        System.out.println("*** zetaMidMeanOdd " + zetaMidMean[1]/N);
        System.out.println("*** zetaGramMeanOdd " + zetaGramMean[1]/N);
        System.out.println("*** zetaMidMeanEven " + zetaMidMean[0]/N);
        System.out.println("*** zetaGramMeanEven " + zetaGramMean[0]/N);
        
        Interpolate.consolidatedF();
        
        File file = new File("out/gSeriesE12/fseries.csv");
        BufferedReader calcInput = new BufferedReader(new FileReader(file));
        calcInput.readLine();
        File outfile = new File("out/gSeries" + Interpolate.prefix
        		+ "/fseriesEstimate.csv");
        PrintWriter out = new PrintWriter(outfile);
        out.println("n-3945951431270L,  f.real, f.im, zeta, diff");
        boolean doCompare = Interpolate.prefix.equals("E12");
        for (int i = 1; i <= N; i++) {
            String diff = ", , ";
        	//Gram
        	if(doCompare) {
				String cf = calcInput.readLine();
				String[] parsed = cf.split("[,\\s]+");
				double cfImF = Double.parseDouble(parsed[2]);
				if (i > 1) {
					diff += nf.format(Math.abs(cfImF - Interpolate.fAtBeta[i - 1][1]));
				}
        	}
            out.println((i+2) + ", " + nf.format(Interpolate.fAtBeta[i-1][0])
                + ", " + nf.format(Interpolate.fAtBeta[i-1][1]) + 
                ", " + nf.format(zetaGram[i-1]) + diff
                );
            //mid
        	if(doCompare) {
				String cf = calcInput.readLine();
				String[] parsed = cf.split("[,\\s]+");
	            double cfReF = Double.parseDouble(parsed[1]);
	            diff = ", ";
	            if(i<N) {
	            	diff += nf.format(Math.abs(cfReF-Interpolate.imFmid[i-1][0]))
	            			+ ",";
	            }
        	}
            out.println( ", " + nf.format(Interpolate.imFmid[i-1][0])
                + ", " + nf.format(Interpolate.imFmid[i-1][1]) 
                + ", " + nf.format(zetaMid[i-1]) + diff);
        }
        out.close();
        calcInput.close();
    }
    
    @Test @Ignore
    public void testConvergence() {
        final double z0 = 251.306919355202, z1 = 251.62429374748942;
        final double d0 = -51.032529476335846, d1 = 174.60414605716676;
        final double max = -11.500212986748618;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        double gram = 251.3291305486841;
        incr = (gram-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        assertEquals(-1.283921548, poly4.eval(gram), 0.1);
        System.out.println();
        gram = 251.57272959458808;
        incr = (gram-poly4.positionMax)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = poly4.positionMax + i*incr;
            System.out.println("n 33, " + nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        assertEquals(-7.45540885886704, poly4.eval(gram), 0.2);
    }

    @Test @Ignore
    public void testMain() {
        final double z0 = 251.62429374748942, z1 = 252.0346709114582;
        final double d0 = 174.60414605716676, d1 = -207.33685149858513;
        final double max = 30.964994982487042;
        //244.2006260 zetaFromRiemann -0.7453399242908442
        //244.2627835 zetaFromRiemann -1.2321436376486554
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        System.out.println();
        System.out.println(nf.format(poly4.positionMax) 
                + ", " + nf.format(poly4.der(poly4.positionMax))
                + ", " + nf.format(poly4.eval(poly4.positionMax))
        );
        double gram = 251.8163286404921;
        incr = (gram-251.62429374748942)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = 251.62429374748942 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        assertEquals(30.597, poly4.eval(gram), 0.1);
    }

}
