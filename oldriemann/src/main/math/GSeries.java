package math;

import java.math.BigDecimal;
import java.util.Arrays;

import javafx.util.Pair;
import riemann.Gram;
import riemann.Interpolate;
import riemann.Riemann;

/**
 * Calculates the G-series useful in the study of Riemann zeta
 * function.
 * zeta(1/2 + it) = exp(−i*theta(t))Z(t)
 * tau = sqrt(t/(2*pi))
 * Z(t) = Real(exp(−i*theta(t))F(1,floor(tau); t)) + R(t)
 *
 * @author oshanker
 *
 */
public class GSeries {
    final int k0;
    public final int k1;
    /**
     * rotation from F to G:
     * G(t) = exp(−i*alpha*t)F(t)
     */
    final double alpha;

    /**
     * Store g at n*beta (k is k0 to k1)
     */
    public double[][] gAtBeta;
    private final double beta;
    public final double spacing;
    private final double gamma;
    private double argalphaBase;
    public final double begin;
    public final double  dsqrtArg1;
    double tbase;
    public double basesqrtArg1;
    double lnsqrtArg1;
    double basetheta;

    /**
     * central index into gAtBeta for t0
     */
    public int midIdx;

    /**
     * create GSeries with pre-calculated gAtBeta
     * @param k0
     * @param k1
     * @param offset
     * @param begin
     * @param incr
     * @param gAtBeta
     */
    public GSeries(int k0, int k1, BigDecimal offset, double begin, double incr, double[][] gAtBeta){
        this(k0, k1, offset, begin, incr);
        this.gAtBeta = gAtBeta;
    }

    /**
     * Calculate GSeries at R points
     * @param k0
     * @param k1
     * @param offset
     * @param begin
     * @param incr
     * @param R
     */
    public GSeries(int k0, int k1, BigDecimal offset, double begin, double incr, int R){
        this(k0, k1, offset, begin, incr);
        BigDecimal tBaseBD = new BigDecimal(begin, Gram.mc).add(
                offset, Gram.mc);
        gAtBeta = calculateGSeries( R, tBaseBD);
    }

    /**
     * initialize basic constants
     * @param k0
     * @param k1
     * @param offset
     * @param begin
     * @param incr
     */
    public  GSeries(int k0, int k1, BigDecimal offset, double begin, double incr) {
        this.k0 = k0;
        this.begin = begin;
        this.spacing = incr;
        beta = Math.PI/spacing;
        BigDecimal tBaseBD = new BigDecimal(begin, Gram.mc).add(
                offset, Gram.mc);
        /////
        dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tBaseBD.doubleValue()));
        tbase = tBaseBD.doubleValue();
        BigDecimal t2 = tBaseBD.divide(Gram.bdTWO);
        BigDecimal sqrtArg1 = Gram.sqrt(tBaseBD.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
        basesqrtArg1 = sqrtArg1.doubleValue();
        this.k1 = k1==0?(int)basesqrtArg1:k1;

        BigDecimal alphaBD = (Gram.log(this.k0).add(Gram.log(this.k1))).divide(Gram.bdTWO, Gram.mc);
        alpha = alphaBD.doubleValue();
        argalphaBase = tBaseBD.multiply(alphaBD, Gram.mc).remainder(Gram.pi_2).doubleValue();

        BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, Gram.mc);
        lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
        basetheta = tBaseBD.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
                .subtract(Gram.pi8, Gram.mc).remainder(Gram.pi_2).doubleValue();
        double tau = (Math.log(this.k1) - Math.log(this.k0))/2.0;
        gamma = beta -tau;
    }

    /**
     * initialize the g-series for the given range of terms.
     * precalculates the 1/sqrt(k) and ln(k) terms.
     * @param k0
     * @param k1
     */
    GSeries(int k0, int k1, int n0, int n1) {
        this.k0 = k0;
        this.k1 = k1;
        int R = n1-n0+1;
        gAtBeta = new double[R][2];
        alpha = (Math.log(k0)+ Math.log(k1))/2.0;
        double tau = (Math.log(k1) - Math.log(k0))/2.0;
        double lambda = 2.0d;
        beta = lambda*tau;
        spacing = Math.PI/beta;
        this.begin = spacing*n0;
        gamma = beta -tau;
        dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*n0*spacing));
        for (int i = 0; i < n1-n0+1; i++) {
            double t = (i+n0)*spacing;
            gAtBeta[i] = gSeriesForSmallT(t);
        }
    }

    public void setgAtBeta(double[][] gAtBeta) {
        this.gAtBeta = gAtBeta;
    }

    /**
     Increment gAtBeta at given index by input value
     */
    public void incrementGValueAtIndex(int index, double[] value) {
        this.gAtBeta[index][0] += value[0];
        this.gAtBeta[index][1] += value[1];
    }

    public void incrementGValuesAtIndices(int startIndex, double[][] valuesToIncrement) {
        for (int i = 0; i < valuesToIncrement.length; i++) {
            this.gAtBeta[i+startIndex][0] += valuesToIncrement[i][0];
            this.gAtBeta[i+startIndex][1] += valuesToIncrement[i][1];
        }
    }
    
    public void decrementGValuesAtIndices(int startIndex, double[][] valuesToIncrement) {
        for (int i = 0; i < valuesToIncrement.length; i++) {
            this.gAtBeta[i+startIndex][0] -= valuesToIncrement[i][0];
            this.gAtBeta[i+startIndex][1] -= valuesToIncrement[i][1];
        }
    }
    
    public void incrementGValuesAtIndices(int startIndex, double[] valuesToIncrement) {
        for (int i = 0; i < valuesToIncrement.length/2; i++) {
            this.gAtBeta[i+startIndex][0] += valuesToIncrement[2*i];
            this.gAtBeta[i+startIndex][1] += valuesToIncrement[2*i+1];
        }
    }
    
    public void decrementGValuesAtIndices(int startIndex, double[] valuesToIncrement) {
        for (int i = 0; i < valuesToIncrement.length/2; i++) {
            this.gAtBeta[i+startIndex][0] -= valuesToIncrement[2*i];
            this.gAtBeta[i+startIndex][1] -= valuesToIncrement[2*i+1];
        }
    }

    public  double evaluateDer(double zero, final int initialPadding) {
        double delta = 0.001*spacing;
        return evalDer(zero, initialPadding, delta);
    }

    public double evalDer(double zero, int initialPadding, double delta) {
        double zetaplus = evaluateZeta(zero + delta, initialPadding);
        double zetaminus = evaluateZeta(zero - delta, initialPadding);
        return (zetaplus-zetaminus)/(2* delta);
    }

    public double[] doubleDer(double zero, int initialPadding, double mid, double delta) {
        double zetaplus = evaluateZeta(zero + delta, initialPadding);
        double zetaminus = evaluateZeta(zero - delta, initialPadding);
        double doubleder = (zetaplus + zetaminus - 2 * mid) / (delta * delta);
        double der =  (zetaplus-zetaminus)/(2* delta);
        return new double[]{der, doubleder};
    }

    public  double evaluateZeta(double zero, final int initialPadding) {
        double[] gFromBLFI = diagnosticBLFISumWithOffset( zero, 4,
                initialPadding, 1.6E-9, false);
        double zeta = riemannZeta(gFromBLFI, zero);
        return zeta;
    }


    public double riemannZeta(double[] g, double tincr){
        return fAndZ(g, tincr).getValue();
    }

    /**
     * Derivative
     * @param zero
     * @param initialPadding
     * @return
     */
    public Pair<double[], Double> der(double zero, final int initialPadding) {
        double delta = 0.001 * spacing;
        double[] gFromBLFI1 = diagnosticBLFISumWithOffset(zero + delta, 4, initialPadding, 1.6E-9, false);
        Pair<double[], Double> fAndZplus = fAndZ(gFromBLFI1, zero + delta);
        double zetaplus = fAndZplus.getValue();
        gFromBLFI1 = diagnosticBLFISumWithOffset(zero - delta, 4, initialPadding, 1.6E-9, false);
        Pair<double[], Double> fAndZminus = fAndZ(gFromBLFI1, zero - delta);
        double zetaminus = fAndZminus.getValue();
        double der = (zetaplus - zetaminus) / (2 * delta);
        double fr = (fAndZplus.getKey()[0] - fAndZminus.getKey()[0]) / (2 * delta);
        double fi = (fAndZplus.getKey()[1] - fAndZminus.getKey()[1]) / (2 * delta);
        double[] fprime = {fr, fi};
        return new Pair<double[], Double>(fprime, der);
    }

    public Pair<double[], Double> fAndZ(double[] g, double tincr){
        double theta = theta(tincr);
        tincr -= begin;
        double predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
        double[] fAtBeta = new double[2];
        double argalphat = (argalphaBase + alpha*tincr)%(2*Math.PI);
        double cos = Math.cos(argalphat);
        double sin = Math.sin(argalphat);
        //calculate f from g
        fAtBeta[0] = cos*g[0] - sin*g[1];
        fAtBeta[1] = sin*g[0] + cos*g[1];
        double rotatedSum = 2*( Math.cos(theta)*fAtBeta [0]+Math.sin(theta)*fAtBeta[1]);
        double correction = GSeries.correction( predictedSqrtArg1);
        double zeta = rotatedSum + correction;
        return new Pair<double[], Double>(fAtBeta, zeta);
    }

    public double theta(double tincr) {
        tincr -= begin;
        double theta = (basetheta + lnsqrtArg1*tincr
                +tincr*tincr/(4*tbase))%(2*Math.PI);
        return theta;
    }

    public  double correctionAtT( double tincr ) {
        tincr -= begin;
        double predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
        return GSeries.correction( predictedSqrtArg1);
    }

    public static double correction( double sqrtArg1 ) {
        long N = (long)sqrtArg1;
        double p = sqrtArg1-N;
        double fourthRoot = Math.sqrt(sqrtArg1);
        double c0formula = Math.cos(2 * Math.PI * (p * p - p - 1.0 / 16))
                / (fourthRoot * Math.cos(2 * Math.PI * p));
        double c1 = Riemann.c1coeff(fourthRoot, p);
        double R = c0formula + c1;
        if (N % 2 == 0) {
            R = -R;
        }
        return R;
    }


    /**
     * Evaluate gSeries without calling cos and sin functions more often than
     * necessary.
     * Offset is placeholder, it will be extended to BigDouble for large height
     * calculations.
     */
    private final  double[][] calculateGSeries( int R, BigDecimal tBase){
        double[][] gAtBeta = fSeries(k0, k1, spacing, R, tBase);
        //now ratate the f to g.
        rotateFtoG(gAtBeta);

        return gAtBeta;
    }

    public void rotateFtoG( double[][] gAtBeta) {
        int R  = gAtBeta.length;
        double costalpha = Math.cos(argalphaBase);
        double sintalpha = Math.sin(argalphaBase);
        double cosdalpha = Math.cos(spacing*alpha);
        double sindalpha = Math.sin(spacing*alpha);
        for (int j = 0; j < R; j++) {
            double tmp = gAtBeta[j][0]*costalpha + gAtBeta[j][1]*sintalpha;
            gAtBeta[j][1] = -gAtBeta[j][0]*sintalpha + gAtBeta[j][1]*costalpha;
            gAtBeta[j][0] = tmp;
            //now set values for next t
            double tmpCos = costalpha*cosdalpha - sintalpha*sindalpha;
            sintalpha = sintalpha*cosdalpha + costalpha*sindalpha;
            costalpha = tmpCos;
        }
    }

    public static double[][] fSeries(int k0, long k1, double incr, int R, BigDecimal tBase) {
        double[][] fAtBeta = new double[R][2];
        for (int i = k0; i <= k1; i++) {
            //evaluate one term in the series, for all t.
            double coeff = 1/Math.sqrt(i);;
            double argi = tBase.multiply(Gram.log(i), Gram.mc).remainder(Gram.pi_2).doubleValue();
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

    /**
     * Calculate the gSeries for arg t.
     * @param t
     * @return
     */
    double[] gSeriesForSmallT(double t){
        double[] g = new double[2];
        double f0 = 0, f1 = 0;
        for (int i = k0; i <= k1; i++) {
            f0 += Math.cos(t*Math.log(i))/Math.sqrt(i);
            f1 += Math.sin(t*Math.log(i))/Math.sqrt(i);
        }
        g[0] = Math.cos(alpha*t)*f0 + Math.sin(alpha*t)*f1;
        g[1] = Math.cos(alpha*t)*f1 - Math.sin(alpha*t)*f0;
        return g;
    }

    private  double[] testblfiSumSmallT( double t0, int M) {
        double[] directEval = gSeriesForSmallT(t0);
        double[] sum = new double[]{0,0};
        //int midIdx = (int) (t0/spacing);
        int midIdx = (int) ((t0-begin)/spacing);
        SUM:
        for (int term = 0; term < 100; term++) {
            for (int j = 0; j < 2; j++) {
                int i = midIdx + (j==0?-term:(term+1));
                if(i>=gAtBeta.length || i < 0){
                    System.out.println("breaking " + i);
                    break SUM;}
                double t = begin + (i)*spacing;
                double harg = gamma*(t0-t)/M ;
                double h = Math.pow( Math.sin(harg)/harg, M);
                double sarg = beta*(t0-t) ;
                double sin = Math.sin(sarg)/sarg;
                sum[0] += gAtBeta[i][0]*h*sin;
                sum[1] += gAtBeta[i][1]*h*sin;
            }
            System.out.println((term) + " : " + Arrays.toString(sum)
                    + " : " + (Math.abs(sum[0] - directEval[0]) + " : " + Math.abs(sum[1] - directEval[1])));
        }
        return sum;
    }

    /**
     * Find the multiplication factor at index i,
     * for arg value t0
     * @param i
     * @param t0
     * @param M power
     * @return
     */
    double factorAtIndex(int i, double t0, int M ){
        double t = begin + (i )*spacing;
        double harg = gamma*(t0-t)/M ;
        if(harg == 0.0d){
            return 1.0;
        }
        double h = Math.pow( Math.sin(harg)/harg, M);
        double sarg = beta*(t0-t) ;
        double sin = Math.sin(sarg)/sarg;
        if(!Double.isFinite(sin)){
            System.out.println(sarg + " i " + i);
            throw new IllegalArgumentException("sin NaN");
        }
        return h*sin;
    }

    public  double[] diagnosticBLFISumWithOffset( double t0, int M,
                                                  int terms, double tolerance, boolean print) {
        double[] sum = new double[]{0,0};
        double[] oldSum = new double[]{0,0};
        // central index into gAtBeta for t0
        midIdx = (int) ((t0-begin)/spacing);
        SUM:
        for (int term = 0; term < terms; term++) {
            oldSum[0] = sum[0];
            oldSum[1] = sum[1];
            for (int j = 0; j < 2; j++) {
                int i = midIdx + (j==0?-term:(term+1));
                if(i>=gAtBeta.length || i < 0){
                    throw new IllegalStateException(t0 + ", index out of bounds " + i);
                }
                double t = begin + (i)*spacing;
                double harg = gamma*(t0-t)/M ;
                double h = Math.pow( Math.sin(harg)/harg, M);
                double sarg = beta*(t0-t) ;
                if(sarg == 0.0d){
                    sum[0] = gAtBeta[i][0];
                    sum[1] = gAtBeta[i][1];
                    return sum;
                }
                double sin = Math.sin(sarg)/sarg;
                if(!Double.isFinite(sin)){
                    System.out.println(sarg + " i " + i);
                    throw new IllegalArgumentException("sin NaN");
                }
                sum[0] += gAtBeta[i][0]*h*sin;
                sum[1] += gAtBeta[i][1]*h*sin;
            }
            double change = Math.abs(sum[0] - oldSum[0]) +  Math.abs(sum[1] - oldSum[1]);
            if(print && (term%5==0 ||(change < tolerance))){
                System.out.println((term) + " : " + Arrays.toString(sum)
                        + " : " + (Math.abs(sum[0] - oldSum[0]) + " : " + Math.abs(sum[1] - oldSum[1])));
            }
            if(change < tolerance){break;}
        }
        return sum;
    }

    double[] blfiSumWithOffsetSmallT( double t0, int M) {
        double[] sum = new double[]{0,0};
        midIdx = (int) ((t0-begin)/spacing);
        SUM:
        for (int term = 0; term < 8; term++) {
            for (int j = 0; j < 2; j++) {
                int i = midIdx + (j==0?-term:(term+1));
                if(i>=gAtBeta.length || i < 0){break SUM;}
                double t = begin + (i)*spacing;
                double harg = gamma*(t0-t)/M ;
                double h = Math.pow( Math.sin(harg)/harg, M);
                double sarg = beta*(t0-t) ;
                double sin = Math.sin(sarg)/sarg;
                sum[0] += gAtBeta[i][0]*h*sin;
                sum[1] += gAtBeta[i][1]*h*sin;
            }
        }
        return sum;
    }

    public static void main(String[] args){
        int k0 = 10, k1=100;
        int N = 25;
        int minIndex = 5;
        GSeries x = new GSeries(k0, k1,minIndex, minIndex+N-1);
        System.out.println("pi/beta " + x.spacing);
        double[] offsets = { 0.3, 0.5, 0.7};
        for (int i = 0; i < offsets.length; i++) {
            double t0 = (minIndex+N/2+offsets[i])*x.spacing;
            double[] sum = x.testblfiSumSmallT( t0, 3);
            System.out.println(t0 + " sum " + Arrays.toString(sum) + ": " + Arrays.toString(x.gSeriesForSmallT(t0)));
        }
    }

}
