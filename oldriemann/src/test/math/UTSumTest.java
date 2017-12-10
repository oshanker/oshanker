package math;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.Arrays;

import org.junit.Test;

import math.UTSum.A1A0r;
import riemann.Gram;

public class UTSumTest {
    public enum SOURCE{GSERIES, UTSUM};
    @Test
    public void testZeroLargeOffset() {
        long offset = (long) 1.0E12;
        double[] begin = {243.8749480149, 1436233.106281030331450810};
//        evaluateZeta(offset, begin, SOURCE.GSERIES);
        evaluateZeta(offset, begin, SOURCE.UTSUM);
    }

    @Test
    public void testTlni() {
        long currentK = 1L<<13;
        BigDecimal tBase = BigDecimal.valueOf(1L<<42).divide(Gram.log(2),  UTSum.mc);
        BigDecimal tlni = tBase.multiply(Gram.log(currentK),  UTSum.mc);
        //57174604644352
        System.out.println(tBase + " multiply " + Gram.log(currentK) + ", " + tlni);
        tlni = tlni.subtract(new BigDecimal(tlni.toBigInteger()));
        System.out.println(".." + tlni);
        double argi = tlni.doubleValue();
        
        long h = (1l<<3);
        //57180797948788.24715053910002
        BigDecimal tlnih = tBase.multiply(Gram.log(currentK+h), UTSum.mc);
        System.out.println( "Gram.log(K+h) " + Gram.log(currentK+h) + ", " + tlnih);
        tlnih = tlnih.subtract(new BigDecimal(tlnih.toBigInteger()));
        System.out.println(".." + tlnih);
        double argih = tlnih.doubleValue();
        
        System.out.println(argi+ ", " + argih +  ", h " + h + ", " +  (argih-argi));
        A1A0r[] coeff = evaluateCoefficients(currentK, tBase, h);

        double sum = performSum(coeff, currentK,  h, tBase.doubleValue());
        System.out.println( "sum " + sum);
        assertEquals(argih, argi + sum, 1.0E-14);
    }

    public double performSum(A1A0r[] coeff, long currentK, long h, double tBaseDbl) {
        double sum = 0;
        
        double uKincr1 = UTSum.calculateIncr1(coeff[0],   h);
        sum += uKincr1;
        System.out.println("** uKincr " + uKincr1  );

        double uKincr2 = UTSum.calculateIncr2(coeff[1],   h, 2);
        sum += uKincr2;
        System.out.println("** uKincr " + uKincr2 + ", " +  (sum)  );

        double uKincr = UTSum.calculateIncr2(coeff[2],   h, 3);
        sum += uKincr;
        System.out.println("** uKincr " + uKincr + ", " +  (sum) );
        
        double rho = ((double)h)/currentK;
        double term = -tBaseDbl*Math.pow(rho, 4)/(4);
        double remaining = term;
        for(int j = 5; j <12; j++){
            term = -rho*term*(j-1)/j;
            remaining += term;
            if(Math.abs(term)<1.0E-15){
                break;
            }
        }
        
        //4: 1.442695040888963, 5: 0.001127105500695, 6: 0.000000917240805
        sum += remaining;
        return sum;
    }

    public A1A0r[] evaluateCoefficients(long currentK, BigDecimal tBase, long h) {
        A1A0r[] coeff = new A1A0r[3];
        
        BigDecimal UK = tBase.divide(BigDecimal.valueOf(currentK));
        //6196328018.719480601951602
        System.out.println( UK + " UK*h " + UK.multiply(BigDecimal.valueOf(h)) );
        coeff[0] = UTSum.evalA1A0(UK);

        UK = UK.divide(BigDecimal.valueOf(2*currentK));
        BigDecimal ukh2 = UK.multiply(BigDecimal.valueOf(h*h));
        //3025550.790390371387672
        System.out.println( UK + ", "  + " UK2*h*h " + ukh2 );
        coeff[1] = UTSum.evalA1A0(UK);
        
        UK = UK.multiply(Gram.bdTWO).divide(BigDecimal.valueOf(3*currentK), UTSum.mc);
        //1969.759629160398039
        System.out.println(  UK + ", "  + " UK*h*h*h " 
        + UK.multiply(BigDecimal.valueOf(h*h*h)) );
        coeff[2] = UTSum.evalA1A0(UK);
        return coeff;
    }
    

    public void evaluateZeta(long offset, double[] begin, SOURCE source) {
        int k0 = 1, k1=0;
        int R = 2;
        double lnsqrtArg1 = 0;
        double basetheta = 0;
        double dsqrtArg1 = 0;
        double tbase = 0;
        double basesqrtArg1 = 0;
        double[][] fAtBeta = null;
        for (int i = 0; i < begin.length; i++) {
            double tincr =  (begin[i]-begin[0]) ; 
            BigDecimal tval = new BigDecimal(begin[i], Gram.mc).add(
                    BigDecimal.valueOf(offset), Gram.mc);
            double predictedSqrtArg1 = 0;
            double theta = 0;
            if(i == 0 ){
                dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
                tbase = tval.doubleValue();
                BigDecimal t2 = tval.divide(Gram.bdTWO);
                BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
                basesqrtArg1 = sqrtArg1.doubleValue();
                k1 = (int)basesqrtArg1;
                BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, Gram.mc);
                lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
                
                long init= System.currentTimeMillis();
                switch (source) {
                case GSERIES:
                    fAtBeta = GSeries.fSeries(k0, k1, begin[1]-begin[0], R, tval);
                    break;

                default:
                    fAtBeta = UTSum.fSeries(k0, k1, begin[1]-begin[0], R, tval);
                    break;
                }
                long end = System.currentTimeMillis();
                System.out.println("calc for " + R + ": " + (end - init) + "ms");
                
                theta = tval.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
                        .subtract(Gram.pi8, Gram.mc).remainder(Gram.pi_2).doubleValue();
                basetheta = theta;
                predictedSqrtArg1 = basesqrtArg1 ;
            } else {
                theta = (basetheta + lnsqrtArg1*tincr
                        +tincr*tincr/(4*tbase))%(2*Math.PI);
                predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
            }
            double rotatedSum = 2*( Math.cos(theta)*fAtBeta[i][0]+Math.sin(theta)*fAtBeta[i][1]);
            double correction = GSeries.correction( predictedSqrtArg1);
            double zeta = rotatedSum + correction;
            System.out.println("f  : " + Arrays.toString(fAtBeta[i])
               + " theta " + theta + "\n rotatedSum " + rotatedSum
               + " zeta " + zeta);
            System.out.println("sqrtArg1[i].doubleValue() " + predictedSqrtArg1 + " correction " + correction );
            assertTrue(i + " ", Math.abs(zeta)  < 5.0E-7);
        }
    }


    @Test
    public void testCalculateIncr2() {
        long h = (1 << 11);
        assertEquals(UTSum.calculateIncr2(
                UTSum.evalA1A0(BigDecimal.valueOf((1l << 32) + (1L << 20) + 0.5d).divide(UTSum.bBD2, UTSum.mc)),
                (1L << 5), 2), -(1.0d + (1L << 21) + (1L << 33)) / (1L << 54), 1.0E-22);

        BigDecimal UK = BigDecimal.valueOf((1l << 32) + (1L << 30) + 0.5d).divide(UTSum.bBD2, UTSum.mc);

        System.out.println("\n" + UK + " UK*h " + UK.multiply(BigDecimal.valueOf(h)));
        A1A0r coeff1 = UTSum.evalA1A0(UK);
        double uKincr1 = UTSum.calculateIncr1(coeff1, h);
        System.out.println("** uKincr " + uKincr1 + ", ");
        System.out.println("** uKincr " + (1.0d + (1L << 31) + (1L << 33)) / (1L << 53) + "\n");
        assertEquals(uKincr1, (1.0d + (1L << 31) + (1L << 33)) / (1L << 53), 1.0E-22);

        BigDecimal ukh2 = UK.multiply(BigDecimal.valueOf(h * h));
        System.out.println("\n" + UK + ", " + " UK2*h*h " + ukh2);

        double uKincr2 = UTSum.calculateIncr2(coeff1, h, 2);
        System.out.println("** uKincr " + uKincr2 + ", " + (uKincr1 + uKincr2));
        System.out.println("** uKincr " + -(1.0d + (1L << 31) + (1L << 33)) / (1L << 42) + ", ");
        assertEquals(uKincr2, -(1.0d + (1L << 31) + (1L << 33)) / (1L << 42), 1.0E-19);

        System.out.println("\n" + UK + ", " + " UK*h*h*h " + ukh2.multiply(BigDecimal.valueOf(h)));
        double uKincr = UTSum.calculateIncr2(coeff1, h, 3);
        System.out.println("** uKincr " + uKincr + ", " + (uKincr1 + uKincr));
        System.out.println("** ****** " + 1.0d / (1L << 31) + ", ");
        assertEquals(uKincr, 1.0d / (1L << 31), 1.0E-28);
    }

}
