package riemann;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Random;

import org.junit.Ignore;
import org.junit.Test;

import math.GSeries;
import riemann.Riemann.GramInfo;

public class ConjecturesTest {
   //45 works well with -45, but 40 doesn't 
   static MathContext mc = new MathContext(45, RoundingMode.HALF_EVEN);
   static double offset = 1.0E28;
   static BigDecimal gram = new BigDecimal("100.3677724301761067");
   static BigDecimal tval = gram.add(
           BigDecimal.valueOf(offset), mc);
   static double incr = 0.10031507797926817;
   /*
    why do we lose significance?
    26: 0.7935445010 //can exceed long limits
    25: 0.7935445010 cf .7935445018
    24  0.793544471  cf .79354450
    */
   static long b = 1L<<32;
   static long mask32 = b-1;

   static double bdbl  = (double)b;
   static long pow31 = (1L<<31);
   static long mask31 = pow31-1;
   static double pow31dbl = (double)pow31;
   static double bdlSquared = 1.0d + (double)Long.MAX_VALUE;
   //static BigDecimal bBD = BigDecimal.valueOf(b);
   static BigDecimal bBD2 = BigDecimal.valueOf(bdlSquared);
   
   public static double[][] fSeries(int k0, long k1, double incr, int R, BigDecimal tBase) {
        double[][] fAtBeta = new double[R][2];
        tBase = tBase.divide(Gram.pi_2,Gram.mc);
        for (int i = k0; i <= k1; i++) {
            //evaluate one term in the series, for all t.
            double coeff = 1/Math.sqrt(i);;
            BigDecimal tlni = tBase.multiply(Gram.log(i), Gram.mc);
            tlni = tlni.subtract(new BigDecimal(tlni.toBigInteger()));
            double argi = tlni.doubleValue()*2*Math.PI;
            double costlni = Math.cos(argi);
            double sintlni = Math.sin(argi);
            //this speeds up, but do we lose accuracy?
            // no, that is not culprit: R costlni -0.9949659697160186, cf -0.9949659697160186
            double cosdlni = Math.cos(incr*Math.log(i));
            double sindlni = Math.sin(incr*Math.log(i));
            for (int j = 0; j < R; j++) {
                fAtBeta[j][0] += coeff*costlni;
                fAtBeta[j][1] += coeff*sintlni;
                //now set values for next t
                double tmpCos = costlni*cosdlni - sintlni*sindlni;
                sintlni = sintlni*cosdlni + costlni*sindlni;
                costlni = tmpCos;
            }
        }
        return fAtBeta;
    }
    
    
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

                fAtBeta = fSeries(k0, k1, begin[1]-begin[0], R, tval);
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
        BigDecimal tBase = tval.divide(Gram.pi_2,Gram.mc);
        BigDecimal tlni = tBase.multiply(Gram.log(K), Gram.mc);
        System.out.println(tBase + " multiply " + Gram.log(K) + ", " + tlni);
        tlni = tlni.subtract(new BigDecimal(tlni.toBigInteger()));
        System.out.println(".." + tlni);
        double argi = tlni.doubleValue()*2*Math.PI;
        
        long h = (1<<11);
        BigDecimal tlnih = tBase.multiply(Gram.log(K+h), Gram.mc);
        System.out.println(tBase + " multiply " + Gram.log(K+h) + ", " + tlnih);
        tlnih = tlnih.subtract(new BigDecimal(tlnih.toBigInteger()));
        System.out.println(".." + tlnih);
        double argih = tlnih.doubleValue()*2*Math.PI;
        
        System.out.println(argi+ ", " + argih + ", b " + b + ", h " + h);
        
        BigDecimal UK = BigDecimal.valueOf((1l<<32) + (1L<<30) + 0.5d).divide(bBD2,mc);
        
//        BigDecimal UK = tBase.divide(BigDecimal.valueOf(K));
        
        System.out.println( "\n" + UK + " UK*h " + UK.multiply(BigDecimal.valueOf(h)) );
        A1A0r coeff1 = evalA1A0(UK);
        double uKincr1 = calculateIncr1(coeff1,   h);
        System.out.println("** uKincr " + uKincr1 + ", " +  (argih-argi)/(2*Math.PI) );
        System.out.println("** uKincr " + (1.0d + (1L<<31) + (1L<<33))/(1L<<53) + "\n"   );
        assertEquals(uKincr1, (1.0d + (1L<<31) + (1L<<33))/(1L<<53),1.0E-22);
        
        assertEquals(calculateIncr2(evalA1A0(BigDecimal.valueOf((1l<<32) + (1L<<20) + 0.5d).divide(bBD2,mc)), 
                (1L<<5), 2), 
                -(1.0d + (1L<<21) + (1L<<33))/(1L<<54),1.0E-22);
        
//        UK = UK.divide(BigDecimal.valueOf(2*K));
        BigDecimal ukh2 = UK.multiply(BigDecimal.valueOf(h*h));
        System.out.println( "\n" + UK + ", "  + " UK2*h*h " + ukh2 );
        A1A0r coeff2 = evalA1A0(UK);
        int k = 2;
        double uKincr2 = calculateIncr2(coeff2,   h, k);
        System.out.println("** uKincr " + uKincr2 + ", " +  (uKincr1 + uKincr2)  );
        System.out.println("** uKincr " + -(1.0d + (1L<<31) + (1L<<33))/(1L<<42) + ", "   );
        assertEquals(uKincr2, -(1.0d + (1L<<31) + (1L<<33))/(1L<<42),1.0E-19);

        k = 3;
        System.out.println( "\n" + UK + ", "  + " UK*h*h*h " 
        + ukh2.multiply(BigDecimal.valueOf(h)) );
        double uKincr = calculateIncr2(coeff2,   h, k);
        System.out.println("** uKincr " + uKincr + ", " +  (uKincr1 + uKincr) );
        System.out.println("** ****** " + 1.0d/(1L<<31) + ", "  );
        assertEquals(uKincr, 1.0d/(1L<<31),1.0E-28);
    }

    public class A1A0r{
        long A1;
        long A0;
        double r;
    }

    public double calculateIncr1(A1A0r coeff,  long h) {
        
        long A1h = (coeff.A1*h)&mask31;
        double t1 = A1h/pow31dbl;
        
        long A0h = coeff.A0*h;
        //mod of b squared trivial, if there is no overflow.
        double t2 = A0h/bdlSquared;
        
        double t3 = coeff.r*h;
        
        System.out.println(t1 + ", " + t2 + ", " +  t3);
        double uKincr = t1 + t2 + t3;
        return uKincr;
    }

    public double calculateIncr2(A1A0r coeff,  long h, int k) {
        long A1h = (coeff.A1);
        double t3 = (coeff.r);
        long overflow = 0;
        long A0h = coeff.A0;
        
        for(int i = 0; i<k; i++){
            t3 *= h;
            A1h = (A1h*h)&mask31;
            A0h = A0h*h;
            if(A0h>mask32){
                long overflow1 = (((A0h)>>32));
                for(int j = 1; j < (k-i); j++){
                   overflow1 = (overflow1*h)&mask31;
                }
                overflow += overflow1;
                A0h &= mask32; //max 2^32-1
            }
        }

        A1h = (A1h + overflow)&mask31;
        
        double t2 = A0h/bdlSquared;
        double t1 = A1h/pow31dbl;
        
        System.out.println(t1 + ", " + t2 + ", " +  t3);
        double uKincr = (t1 + t2 + t3);
        if(k%2==0){uKincr = -uKincr;}
        return uKincr;
    }


    public A1A0r evalA1A0(BigDecimal UK) {
        A1A0r coeff = new A1A0r();
        BigDecimal UKfrac = UK.subtract(new BigDecimal(UK.toBigInteger()));
        BigDecimal UKfracNorm = UKfrac.multiply(bBD2);
        BigDecimal intValue = new BigDecimal(UKfracNorm.toBigInteger());
        long UKNormint = intValue.longValue();

        coeff.r = UKfracNorm.subtract(intValue).divide(bBD2,mc).doubleValue();
        coeff.A0 = UKNormint&mask32;
        coeff.A1 = (UKNormint-coeff.A0)>>32;
        
        System.out.println(" A1 " +  coeff.A1 + ", A0 " + coeff.A0  + ", r " + coeff.r);
        System.out.println(" test " + (coeff.A1/pow31dbl +  (coeff.A0)/bdlSquared  +  coeff.r));
        return coeff;
    }

    @Test
    public void testScaleFactors() {
        
        System.out.println("gram " + gram);
        checkGramValue( tval,  incr);
    }
    
    void checkGramValue(BigDecimal tval, double incr) {
//        incr = 0.1921410553288139095;
//        BigDecimal tval = new BigDecimal(192.2043554309546, mc).add(
//                BigDecimal.valueOf(1.0E15), mc);
        BigDecimal xx = tval.divide(Gram.pi_2, mc);
        //39894228040143.2677939946061937820389927342029 65420163328854914850403286571752375 
        //39894228040143.2677939946061937820389927342029 24157013297165575084867900724971497
        BigDecimal sqrtArg1 = sqrt(xx, mc, 1.0E-45);
        BigDecimal fourthrootArg1 = sqrt(sqrtArg1, mc, 1.0E-45);
        int k1 = fourthrootArg1.intValue();
        Gram.initLogVals(k1);
        BigDecimal thetaPi = theta(tval, fourthrootArg1, mc);
        System.out.println(thetaPi);
        
        tval = tval.add(new BigDecimal(Double.toString(incr)));
        sqrtArg1 = sqrt(tval.divide(Gram.pi_2, mc), mc, 1.0E-45);
        fourthrootArg1 = sqrt(sqrtArg1, mc, 1.0E-45);
        BigDecimal thetaPi2 = theta(tval, fourthrootArg1, mc);
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
