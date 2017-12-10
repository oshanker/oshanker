package riemann;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import math.GSeries;
import math.UTSum;
import math.UTSum.A1A0r;

public class ConjecturesTest {
   //45 works well with -45, but 40 doesn't 
   static double offset = 1.0E28;
   static BigDecimal gram = new BigDecimal("100.3677724301761067");
   static BigDecimal tval = gram.add(
           BigDecimal.valueOf(offset), UTSum.mc);
   static double incr = 0.10031507797926817;
   //@Test //needs fixing (large sums)
    public void testZeroLargeOffset() {
        MathContext mc = new MathContext(50, RoundingMode.HALF_EVEN);
        double[][] fAtBeta = null;
        double[] begin = {100.437512887104287873, 100.464843234223048518};
        int k0 = 1;
        long k1;
        int R = 2;
        double lnsqrtArg1 = 0;
        double basetheta = 0;
        double dsqrtArg1 = 0;
        double basesqrtArg1 = 0;
        BigDecimal offset =  BigDecimal.valueOf(1.0E28);
        for (int i = 0; i < begin.length; i++) {
            double tincr =  (begin[i]-begin[0]) ; 
            BigDecimal tval = new BigDecimal(begin[i], mc).add(
                    offset, mc);
            double predictedSqrtArg1 = 0;
            double theta = 0;
            if(i == 0 ){
                dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
                BigDecimal xx = tval.divide(Gram.pi_2, mc);
                BigDecimal sqrtArg1 = sqrt(xx, mc, 1.0E-38);
                BigDecimal fourthrootArg1 = sqrt(sqrtArg1, mc, 1.0E-38);
                Gram.initLogVals((int)fourthrootArg1.intValue()/2);
                k1 = (sqrtArg1.longValue());
                System.out.println(sqrtArg1);
                System.out.println("k1 " + k1);
                k1 = (long) 3989420;
                // a billion seconds
                BigDecimal thetaPi = theta(tval, fourthrootArg1, mc);
                System.out.println(thetaPi);

                BigDecimal t2 = tval.divide(Gram.bdTWO);
                basesqrtArg1 = sqrtArg1.doubleValue();
                BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, mc);
                lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
                long init= System.currentTimeMillis();

                fAtBeta = UTSum.fSeries(k0, k1, begin[1]-begin[0], R, tval);
                long end = System.currentTimeMillis();
                System.out.println("evaluateWithOffset calc for " + k1 + ": " + (end - init) + "ms");
                //theta should be 2.819633653651107
                theta = tval.multiply(lnsqrtArg1BD, mc).subtract(t2, mc)
                        .subtract(Gram.pi8, mc).remainder(Gram.pi_2).doubleValue();
                basetheta = theta;
                predictedSqrtArg1 = basesqrtArg1 ;
                System.out.println(predictedSqrtArg1);
            } else {
                theta = (basetheta + lnsqrtArg1*tincr)%(2*Math.PI);
                predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
            }
            double rotatedSum = 2*( Math.cos(theta)*fAtBeta[i][0]+Math.sin(theta)*fAtBeta[i][1]);
            double correction = GSeries.correction( predictedSqrtArg1);
            double zeta = rotatedSum + correction;
            System.out.println("f  : " + Arrays.toString(fAtBeta[i])
               + " theta " + theta + " rotatedSum " + rotatedSum
               + " zeta " + zeta);
            System.out.println("sqrtArg1[i].doubleValue() " + predictedSqrtArg1 + " correction " + correction );
            assertTrue(i + " ", Math.abs(zeta)  < 5.0E-7);
        }
    }

    @Test
    public void testTlni() {
        long K = (long) 1.0E8;
        BigDecimal tBase = tval.divide(Gram.pi_2,UTSum.mc);
        BigDecimal tlni = tBase.multiply(Gram.log(K), Gram.mc);
        System.out.println(tBase + " multiply " + Gram.log(K) + ", " + tlni);
        tlni = tlni.subtract(new BigDecimal(tlni.toBigInteger()));
        System.out.println(".." + tlni);
        double argi = tlni.doubleValue()*2*Math.PI;
        
        long h = (1L<<11);
        BigDecimal tlnih = tBase.multiply(Gram.log(K+h), Gram.mc);
        System.out.println(tBase + " multiply " + Gram.log(K+h) + ", " + tlnih);
        tlnih = tlnih.subtract(new BigDecimal(tlnih.toBigInteger()));
        System.out.println(".." + tlnih);
        double argih = tlnih.doubleValue()*2*Math.PI;
        
        System.out.println(argi+ ", " + argih + ", b " + UTSum.b + ", h " + h);
        
        
        BigDecimal UK = tBase.divide(BigDecimal.valueOf(K));
        
        System.out.println( "\n" + UK + " UK*h " + UK.multiply(BigDecimal.valueOf(h)) );
        A1A0r coeff1 = UTSum.evalA1A0(UK);
        double uKincr1 = UTSum.calculateIncr1(coeff1,   h);
        System.out.println("** uKincr " + uKincr1 + ", " +  (argih-argi)/(2*Math.PI) );
        
        
        UK = UK.divide(BigDecimal.valueOf(2*K));
        BigDecimal ukh2 = UK.multiply(BigDecimal.valueOf(h*h));
        System.out.println( "\n" + UK + ", "  + " UK2*h*h " + ukh2 );
        A1A0r coeff2 = UTSum.evalA1A0(UK);

        double uKincr2 = UTSum.calculateIncr2(coeff2,   h, 2);
        System.out.println("** uKincr " + uKincr2 + ", " +  (uKincr1 + uKincr2)  );

        UK = UK.multiply(Gram.bdTWO).divide(BigDecimal.valueOf(3*K));
        System.out.println( "\n" + UK + ", "  + " UK*h*h*h " 
        + UK.multiply(BigDecimal.valueOf(h*h*h)) );
        A1A0r coeff3 = UTSum.evalA1A0(UK);
        double uKincr = UTSum.calculateIncr2(coeff3,   h, 3);
        System.out.println("** uKincr " + uKincr + ", " +  (uKincr1 + uKincr2 + uKincr) );
    }
    
    @Test
    public void testScaleFactors() {
        
        System.out.println("gram " + gram);
        checkGramValue( tval,  incr);
    }
    
    void checkGramValue(BigDecimal tval, double incr) {
//        incr = 0.1921410553288139095;
//        BigDecimal tval = new BigDecimal(192.2043554309546, UTSum.mc).add(
//                BigDecimal.valueOf(1.0E15), UTSum.mc);
        BigDecimal xx = tval.divide(Gram.pi_2, UTSum.mc);
        //39894228040143.2677939946061937820389927342029 65420163328854914850403286571752375 
        //39894228040143.2677939946061937820389927342029 24157013297165575084867900724971497
        BigDecimal sqrtArg1 = sqrt(xx, UTSum.mc, 1.0E-45);
        BigDecimal fourthrootArg1 = sqrt(sqrtArg1, UTSum.mc, 1.0E-45);
        int k1 = fourthrootArg1.intValue();
        Gram.initLogVals(k1);
        BigDecimal thetaPi = theta(tval, fourthrootArg1, UTSum.mc);
        System.out.println(thetaPi);
        
        tval = tval.add(new BigDecimal(Double.toString(incr)));
        sqrtArg1 = sqrt(tval.divide(Gram.pi_2, UTSum.mc), UTSum.mc, 1.0E-45);
        fourthrootArg1 = sqrt(sqrtArg1, UTSum.mc, 1.0E-45);
        BigDecimal thetaPi2 = theta(tval, fourthrootArg1, UTSum.mc);
        System.out.println(thetaPi2);
        
        assertEquals(thetaPi2.subtract(thetaPi).doubleValue(), 1, 1.0E-15);
       /*
baseGram 192.2043554309546 idx 5045354828590534 incr at mid 0.19214105532316908
        double[] begin = {192.309350419702134727, 1921602.793342093316678265};
f  : [1.9524278734805527, 0.2869931460545315] theta 1.716715278515568 rotatedSum 1.1379834185443194E-4 zeta 7.381699365530212E-11
zetaFromRiemann 7.413854952046176E-11
sqrtArg1[i].doubleValue() 1.2615662610102013E7 correction -1.1379826803743828E-4
g  : [-1.23276504691923, 0.32595788950860927] argalphaBase 2.4139476198307377
f  : [0.7037596046854905, -1.063334672364068] theta 0.5846098679289753 rotatedSum 1.1521610619680267E-4 zeta -1.3689409652075171E-8
zetaFromRiemann 4.895887538723708E-9
sqrtArg1[i].doubleValue() 1.2615662622221947E7 correction -1.1522979560645475E-4
         */
    }

    /**
     * This differs from original Gram sqrt in the way the loop values
     * are evaluated. This implementation is more stable than the original implementation, 
     * the original tends to sometimes not terminate.
     * @param x
     * @param mc
     * @param prec
     * @return
     */
    public static BigDecimal sqrt(BigDecimal x, MathContext mc, double prec) {
        double init = Math.sqrt(x.doubleValue());
        BigDecimal next = BigDecimal.valueOf(init);
        BigDecimal diff = null;
        while (true) {
            BigDecimal oldnext = next;
            next = x.divide(next, mc).add(next, mc).divide(Gram.bdTWO,mc);
            diff =  next.subtract(oldnext, mc);
            if(Math.abs(diff.doubleValue()) < prec){ 
                break; 
            }
        }
        return next;
    }
    
    private BigDecimal theta(BigDecimal tval, BigDecimal fourthrootArg1, MathContext mc) {
        BigDecimal lnsqrtArg1BD = Gram.log(fourthrootArg1, mc).multiply(Gram.bdTWO);
//        System.out.println(Gram.pi.divide(lnsqrtArg1BD, mc));
//        BigDecimal thetaPi = tval.multiply((lnsqrtArg1BD.subtract(BigDecimal.ONE.divide(Gram.bdTWO), mc)), mc)
//                .subtract(Gram.pi8, mc).divide(Gram.pi, mc);
        BigDecimal t2 = tval.divide(Gram.bdTWO);
        BigDecimal thetaPi = tval.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
             .subtract(Gram.pi8, Gram.mc).divide(Gram.pi, Gram.mc);
        return thetaPi;
    }

    private static double increment(double t) {
        double incr = Math.PI/Math.log(Math.sqrt(t/(2*Math.PI)));
        System.out.println(incr);
        t += 5000000*incr;
        incr = Math.PI/Math.log(Math.sqrt(t/(2*Math.PI)));
        System.out.println(incr);
        return incr;
    }

	@Test @Ignore
	public void testCalculateDistribution() {
		Random x = new Random(5);
		int size = 1000002;
		char[] series = new char[size];
		int[][] doubletcounts = new int[2][4];
		int[][] tripletcounts = new int[2][8];
		char old1 = '1';
		char old = '1';
		for (int i = 0; i < series.length; i++) {
			float y = x.nextFloat();
			char current = '1';
			if(y<0.2){
				current = i%2==0?'-':'+';
			} else {
				current = i%2==0?'+':'-';
			}
			series[i] = current;
			if(i==0){
				old = current;
				continue;
			}
			int idx = 0;
			switch (old) {
			case '-':
				idx = 0;
				break;
			case '+':
				idx = 2;
				break;

			default:
				throw new IllegalStateException();
			}
			switch (current) {
			case '-':
				idx += 0;
				break;
			case '+':
				idx += 1;
				break;

			default:
				throw new IllegalStateException();
			}
			doubletcounts[(i-1)%2][idx]++;
			if(i==1){
				old1 = old;
				old = current;
				continue;
			}
			switch (old1) {
			case '-':
				idx += 0;
				break;
			case '+':
				idx += 4;
				break;

			default:
				throw new IllegalStateException();
			}
			tripletcounts[(i)%2][idx]++;
			old1 = old;
			old = current;
		}
		//System.out.println(Arrays.toString(series));
		System.out.println(Arrays.toString(Conjectures.descriptions[1]));
		for (int i = 0; i < doubletcounts.length; i++) {
			System.out.println(Arrays.toString(doubletcounts[(i+1)%2]));
		}
		System.out.println(Arrays.toString(Conjectures.descriptions[2]));
		for (int i = 0; i < tripletcounts.length; i++) {
			System.out.println(Arrays.toString(tripletcounts[(i+1)%2]));
		}
	}

}
