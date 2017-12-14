package math;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;

import math.UTSum.SOURCE;
import riemann.Gram;

public class UTSum {
    public static MathContext mc = new MathContext(45, RoundingMode.HALF_EVEN);
    public static final long b = 1L << 32;
    public static final long mask32 = b - 1;

    public static final double bdbl = (double) b;
    public static final long pow31 = (1L << 31);
    public static final long mask31 = pow31 - 1;
    public static final double pow31dbl = (double) pow31;
    public static final double bdlSquared = 1.0d + (double) Long.MAX_VALUE;
    public static final BigDecimal bBD2 = BigDecimal.valueOf(bdlSquared);
    public static final double twoPI = 2*Math.PI;

    public static class A1A0r {
        public long A1;
        public long A0;
        public double r;
        @Override
        public String toString() {
            return "A1A0r [A1=" + A1 + ", A0=" + A0 + ", r=" + r + "]";
        }
    }

    public enum SOURCE{GSERIES, UTSUM, SIGMA1}


    public static double calculateIncr1(A1A0r coeff, long h) {
        long A1h = (coeff.A1 * h) & mask31;
        double t1 = A1h / pow31dbl;

        long A0h = coeff.A0 * h;
        // mod of b squared trivial, if there is no overflow.
        double t2 = A0h / bdlSquared;

        double t3 = coeff.r * h;
        double uKincr = t1 + t2 + t3;
        return uKincr;
    }

    public static double calculateIncr2(final A1A0r coeff, long h, final int k) {
        long A1h = (coeff.A1);
        double t3 = (coeff.r);
        long overflow = 0;
        long A0h = coeff.A0;

        for (int i = 0; i < k; i++) {
            t3 *= h;
            A1h = (A1h * h) & mask31;
            A0h = A0h * h;
            if (A0h > mask32) {
                long overflow1 = (((A0h) >> 32));
                for (int j = 1; j < (k - i); j++) {
                    overflow1 = (overflow1 * h) & mask31;
                }
                overflow += overflow1;
                A0h &= mask32; // max 2^32-1
            }
        }

        A1h = (A1h + overflow) & mask31;

        double t2 = A0h / bdlSquared;
        double t1 = A1h / pow31dbl;

        double uKincr = (t1 + t2 + t3);
        if (k % 2 == 0) {
            uKincr = -uKincr;
        }
        return uKincr;
    }

    public static A1A0r evalA1A0(BigDecimal UK) {
        A1A0r coeff = new A1A0r();
        BigDecimal UKfrac = UK.subtract(new BigDecimal(UK.toBigInteger()));
        BigDecimal UKfracNorm = UKfrac.multiply(bBD2);
        BigDecimal intValue = new BigDecimal(UKfracNorm.toBigInteger());
        long UKNormint = intValue.longValue();

        coeff.r = UKfracNorm.subtract(intValue).divide(bBD2, mc).doubleValue();
        coeff.A0 = UKNormint & mask32;
        coeff.A1 = (UKNormint - coeff.A0) >> 32;

        return coeff;
    }

    public static double[][] fSeries(final double[][] fAtBeta, final int k0, final long k1, final double incr, final int R, 
            BigDecimal tBase) {
        tBase = tBase.divide(Gram.pi_2, mc);
        double tBaseDbl = tBase.doubleValue();
        long K = (long) Math.ceil(Math.pow(tBaseDbl/4.0d, 0.25d));
        System.out.println("K " + K + " k1 " +k1 + " tBaseDbl " + tBaseDbl);

        double uk = 0;
        long i = k0;
        // do loop upto K
        for (; i <= K; i++) {
            //evaluate one term in the series, for all t.
            double coeff1 = 1/Math.sqrt(i);
            
            BigDecimal tlni = tBase.multiply(Gram.log(i),  mc);
            tlni = tlni.subtract(new BigDecimal(tlni.toBigInteger()));
            uk = tlni.doubleValue();
            
            double argi = uk*twoPI;
            
            double costlni = Math.cos(argi);
            double sintlni = Math.sin(argi);

            double cosdlni = Math.cos(incr*Math.log(i));
            double sindlni = Math.sin(incr*Math.log(i));
            
            for (int j = 0; j < R; j++) {
                fAtBeta[j][0] += coeff1*costlni;
                fAtBeta[j][1] += coeff1*sintlni;
                if(j < R - 1){
                    //now set values for next t
                    double tmpCos = costlni*cosdlni - sintlni*sindlni;
                    sintlni = sintlni*cosdlni + costlni*sindlni;
                    costlni = tmpCos;
                }
            }
        }
       
        long updateHmax = 2*K; 
        int h = 0;
        double rho = 0;
        int hmax = 1;
        long currentK = i-1;
        double rhoIncr = 1.0d/currentK;
        A1A0r[] coeff = new A1A0r[3];
        
        BigDecimal UK = tBase.divide(BigDecimal.valueOf(currentK), mc);
        coeff[0] = UTSum.evalA1A0(UK);

        UK = UK.divide(BigDecimal.valueOf(2*currentK), mc);
        coeff[1] = UTSum.evalA1A0(UK);
        
        UK = UK.multiply(Gram.bdTWO).divide(BigDecimal.valueOf(3*currentK), UTSum.mc);
        coeff[2] = UTSum.evalA1A0(UK);

        for (; i <= k1; i++) {
            // There are two special events to check for.
            
            //2. do we need to increment hmax?
            if(i >= updateHmax){
                updateHmax += K;
                hmax += 1;
            }            
            
            //evaluate one term in the series, for all t.
            double coeff1 = 1/Math.sqrt(i);
            double argi_2pi;
            //1. do we need to generate a new uk?
            if(h >= hmax){
                //generate uk
                currentK = i;
                rhoIncr = 1.0d/currentK;
                h = 0;
                rho = 0;
                UK = tBase.divide(BigDecimal.valueOf(currentK), mc);
                coeff[0] = UTSum.evalA1A0(UK);

                UK = UK.divide(BigDecimal.valueOf(2*currentK), mc);
                coeff[1] = UTSum.evalA1A0(UK);
                
                UK = UK.multiply(Gram.bdTWO).divide(BigDecimal.valueOf(3*currentK), UTSum.mc);
                coeff[2] = UTSum.evalA1A0(UK);
                BigDecimal tlni = tBase.multiply(Gram.log(i),  mc);
                tlni = tlni.subtract(new BigDecimal(tlni.toBigInteger()));
                uk = tlni.doubleValue();
                argi_2pi = uk;
            } else {
                //  use expansion and increment h
                h++;
                rho += rhoIncr;
                double sum = 0;
                sum += ((coeff[0].A1 * h) & mask31)/ pow31dbl;

                sum +=  (coeff[0].A0 * h)/ bdlSquared;
                
                if(h > 8000){
                    sum +=  coeff[0].r * h;
                }

               sum += UTSum.calculateIncr2(coeff[1],   h, 2);

                sum += UTSum.calculateIncr2(coeff[2],   h, 3);
                
                double term = -tBaseDbl*Math.pow(rho, 4)/(4);
                double remaining = term;
                for(int j = 5; j <12; j++){
                    term = -rho*term*(j-1)/j;
                    remaining += term;
                    if(Math.abs(term)<1.0E-15){
                        break;
                    }
                }
                
                sum += remaining;
                argi_2pi = uk + sum;
            }
            
            double argi = argi_2pi*twoPI;
            
            double costlni = Math.cos(argi);
            double sintlni = Math.sin(argi);

            double cosdlni = Math.cos(incr*Math.log(i));
            double sindlni = Math.sin(incr*Math.log(i));
            
            for (int j = 0; j < R; j++) {
                fAtBeta[j][0] += coeff1*costlni;
                fAtBeta[j][1] += coeff1*sintlni;
                //now set values for next t
                if(j < R - 1){
                    double tmpCos = costlni*cosdlni - sintlni*sindlni;
                    sintlni = sintlni*cosdlni + costlni*sindlni;
                    costlni = tmpCos;
                }
            }
        }
        System.out.println("hmax " + hmax + ", updateHmax " + updateHmax);
        System.out.println( " currentK " + currentK);
        return fAtBeta;
    }
    
    public static double performSum(A1A0r[] coeff, long currentK, long h, double tBaseDbl) {
        double sum = 0;
        //System.out.println("coeff " + Arrays.toString(coeff) + ",\n h " + h + " currentK " + currentK);
        double uKincr1 = UTSum.calculateIncr1(coeff[0],   h);
        sum += uKincr1;
        //System.out.println("** uKincr " + uKincr1  );

        double uKincr2 = UTSum.calculateIncr2(coeff[1],   h, 2);
        sum += uKincr2;
        //System.out.println("** uKincr " + uKincr2 + ", " +  (sum)  );

        double uKincr = UTSum.calculateIncr2(coeff[2],   h, 3);
        sum += uKincr;
        //System.out.println("** uKincr " + uKincr + ", " +  (sum) );
        
        double rho = ((double)h)/currentK;
        double term = -tBaseDbl*Math.pow(rho, 4)/(4);
        double remaining = term;
        //System.out.println("** term " + term + ", " +  (remaining) );
        for(int j = 5; j <12; j++){
            term = -rho*term*(j-1)/j;
            remaining += term;
            //System.out.println("** term " + term + ", " +  (remaining) );
            if(Math.abs(term)<1.0E-15){
                break;
            }
        }
        
        sum += remaining;
        return sum;
    }

    public static void calculateLargeValue() {
        Gram.initLogVals(100000);
        BigDecimal tFactor = Gram.pi_2.divide(Gram.log(Gram.bdTWO, mc), mc);
        BigDecimal t = tFactor.multiply(BigDecimal.valueOf(3164680085302956L));
        long offset = t.longValue();
        double incr = Math.PI/Math.log(Math.sqrt(offset/(2*Math.PI)));
        double[] begin = {t.subtract(BigDecimal.valueOf(offset)).doubleValue(), 0};
        begin[1] = begin[0]+incr/14.5d;
        double[] zeta = evaluateZeta(offset, begin, UTSum.SOURCE.SIGMA1);
    }

    public static double[] evaluateZeta(long offset, double[] begin, SOURCE source) {
        int k0 = 1, k1=0;
        int R = 2;
        double lnsqrtArg1 = 0;
        double basetheta = 0;
        double dsqrtArg1 = 0;
        double tbase = 0;
        double basesqrtArg1 = 0;
        double[][] fAtBeta = new double[R][6];
        double[] zeta = new double[2];
        for (int i = 0; i < begin.length; i++) {
            double tincr =  (begin[i]-begin[0]) ; 
            BigDecimal tval = new BigDecimal(begin[i], mc).add(
                    BigDecimal.valueOf(offset), mc);
            double predictedSqrtArg1 = 0;
            double theta = 0;
            if(i == 0 ){
                dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
                tbase = tval.doubleValue();
                BigDecimal t2 = tval.divide(Gram.bdTWO);
                BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, mc), mc, 1.0E-21);
                basesqrtArg1 = sqrtArg1.doubleValue();
                k1 = (int)basesqrtArg1;
                BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, mc);
                lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
                
                long init= System.currentTimeMillis();
                switch (source) {
                case GSERIES:
                    fAtBeta = GSeries.fSeries(k0, k1, begin[1]-begin[0], R, tval);
                    break;
    
                case SIGMA1:
                    fAtBeta = fSeries(fAtBeta, k0, k1, begin[1]-begin[0], R, tval);
                    break;
    
               default:
                    fAtBeta = fSeries(fAtBeta, k0, k1, begin[1]-begin[0], R, tval);
                    break;
                }
                long end = System.currentTimeMillis();
                System.out.println("calc for " + R + ": " + (end - init) + "ms");
                
                theta = tval.multiply(lnsqrtArg1BD, mc).subtract(t2, mc)
                        .subtract(Gram.pi8, mc).remainder(Gram.pi_2).doubleValue();
                basetheta = theta;
                predictedSqrtArg1 = basesqrtArg1 ;
            } else {
                theta = (basetheta + lnsqrtArg1*tincr
                        +tincr*tincr/(4*tbase))%(2*Math.PI);
                predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
            }
            double rotatedSum = 2*( Math.cos(theta)*fAtBeta[i][0]+Math.sin(theta)*fAtBeta[i][1]);
            double correction = GSeries.correction( predictedSqrtArg1);
            zeta[i] = rotatedSum + correction;
            System.out.println("f  : " + Arrays.toString(fAtBeta[i])
               + " theta " + theta + "\n rotatedSum " + rotatedSum
               + " *** zeta " + zeta[i]);
            System.out.println("sqrtArg1[i].doubleValue() " + predictedSqrtArg1 + " correction " + correction );
        }
        return zeta;
    }


    public static void main(String[] args) {
        calculateLargeValue();
    }
}
