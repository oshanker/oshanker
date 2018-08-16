package math;

import static org.junit.Assert.*;

import java.math.BigDecimal;

import org.junit.Test;

import math.UTSum.A1A0r;
import riemann.Gram;

public class UTSumTest {
    @Test
    public void testZeroLargeOffset() {
        //1415ms
        long offset = (long) 1.0E12;
        double[] begin = {243.8749480149, 1436233.106281030331450810};
        
        //9870ms
//        long offset = (long) 1.0E15;
//        double[] begin = {192.309350419702134727, 200.446361911804399436};
        
        //27522 
//        long offset = (long) 1.0E16;
//        double[] begin = {179.701837722056461444, 180.110245349184374800};

        //254759ms
        //-ea -Xmx2G -Xms2G 
//        long offset = (long) 1.0E18;
//        double[] begin = {158.997193245871834444, 168.073334879216515686};

        //??
//        long offset = (long) 1.0E20;
//        double[] begin = {142.285612643329601679, 145.088695674076487929};

        Gram.initLogVals(100000);
//        evaluateZeta(offset, begin, SOURCE.GSERIES);
        double[] zeta = UTSum.evaluateZeta(offset, begin, UTSum.SOURCE.UTSUM);
        for (int j = 0; j < zeta.length; j++) {
            assertTrue(j + " ", Math.abs(zeta[j])  < 5.0E-7);
        }

        System.out.println("Gram.logVals.length " + Gram.logVals.length);
    }

    @Test
    public void testLargeValue() {
        Gram.initLogVals(100000);
        long offset = 80587673978410171L;
        double incr = Math.PI/Math.log(Math.sqrt(offset/(2*Math.PI)));
        double[] begin = {0.1558857141254, 0};
        begin[1] = begin[0]+incr;
        double[] zeta = UTSum.evaluateZeta(offset, begin, UTSum.SOURCE.UTSUM);
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

    private double performSum(A1A0r[] coeff, long currentK, long h, double tBaseDbl) {
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

    private A1A0r[] evaluateCoefficients(long currentK, BigDecimal tBase, long h) {
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

}
