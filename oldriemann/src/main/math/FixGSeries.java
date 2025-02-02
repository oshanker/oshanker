package math;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.Arrays;

import riemann.Interpolate;

import static riemann.StaticMethods.*;

public class FixGSeries {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(8);
        nf.setMaximumFractionDigits(8);
        nf.setGroupingUsed(false);
    }

    public static void main(String[] args) throws Exception {
        //oldMain();
        double[] nextValues = new double[]{
              243831.456494008, -22.69554476177354, 1.538114456203189};
        double pointBeingInflunced = nextValues[0];
        int fileIdx = findFile(pointBeingInflunced);
        final int initialPadding = 40;

        double t0 = gramE12[fileIdx][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = getSavedGSeries(t0, offset);
        double[] initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        int midIdxCausingInfluence  = gAtBeta.midIdx;

        System.out.println("initial " + Arrays.toString(initial));
        System.out.println("midIdx " + midIdxCausingInfluence);
        double[] zetaCoeff = changeToZeta(
                gAtBeta, initialPadding, pointBeingInflunced, initial[0], midIdxCausingInfluence, 10);
        System.out.println("zetaCoeff " + Arrays.toString(zetaCoeff));
        double[] derCoeff = changeToDer(gAtBeta, initialPadding, pointBeingInflunced, initial[1],
                midIdxCausingInfluence, 1);
        System.out.println("derCoeff " + Arrays.toString(derCoeff));
        double[][] coefficients = new double[][]{zetaCoeff,derCoeff};
        LinearEquation linearEquation = new LinearEquation(coefficients );
        double[] solution = linearEquation.solve(new double[]{1, 0.5});
        System.out.println("Required g increment " + Arrays.toString(solution));

        gAtBeta.incrementGValueAtIndex(midIdxCausingInfluence, solution);
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        double[] actualIncrement = new double[after.length];
        for (int i = 0; i < actualIncrement.length; i++) {
            actualIncrement[i] = after[i] - initial[i];
        }
        System.out.println("after " + Arrays.toString(after));
        System.out.println("actualIncrement " + Arrays.toString(actualIncrement));
    }

    static double[] evaluateAtT(double t, int initialPadding, GSeries gAtBeta) {
        double zeta = gAtBeta.evaluateZeta(t, initialPadding);
        double der = gAtBeta.evalDer(
              t, initialPadding, 0.00025* gAtBeta.spacing);
        double[] initial = {zeta, der};
        return initial;
    }

    private static void oldMain() throws IOException {
        File file = new File("out/gSeries" + Interpolate.prefix + "/zeros.dat");
        DataInputStream in = Interpolate.dataInputStream( file);
        double[] tmin = new double[4];
        GSeries gSeries = Interpolate.readGSeries();
        System.out.println("gSeries.begin " + gSeries.begin);
        final int initialPadding = 40;
        NumberFormat nf1 = NumberFormat.getInstance();
        nf1.setMinimumFractionDigits(2);
        nf1.setMaximumFractionDigits(2);
        nf1.setGroupingUsed(false);

        int N = 1000022;
        for (int midIdx1 = 39; midIdx1 < 39+N ; midIdx1++)
        {
            for (int i1 = 0; i1 < tmin.length; i1++)
            {
                tmin[i1] = in.readDouble();
            }
            double d = tmin[3];
            double expectedDer = tmin[1];
            double expected = expectedDer<0?-tmin[2]:tmin[2];
            double oldmin = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            double valAtZero = Interpolate.evaluateZeta(tmin[0], initialPadding, gSeries);
            double derAtZero = Interpolate.evaluateDer(tmin[0], initialPadding, gSeries);

            double[] g0incr = FixGSeries.evalGSeriesIncrement(gSeries, midIdx1, initialPadding, tmin);
            gSeries.incrementGValueAtIndex(midIdx1, g0incr);
            double min = Interpolate.evaluateZeta(d, initialPadding, gSeries);

            double deviationMin = min - expected;
            if(Math.abs(deviationMin)>100.0 || midIdx1%100000 == 0){
                System.out.println();
                System.out.println(Arrays.toString(tmin));
                System.out.println("valAtZero " + valAtZero + ", midIdx " + midIdx1
                        + ", derAtZero " + derAtZero);
                System.out.println("eval zero final " + Interpolate.evaluateZeta(tmin[0], initialPadding , gSeries) + ", "
                        + Interpolate.evaluateDer(tmin[0], initialPadding , gSeries)
                        + " cf " + expectedDer);
                System.out.println(nf.format(d) + ", " + nf.format(oldmin) + ", "
                + nf.format(expected) + ", "
                        + nf1.format((oldmin - expected) * 100.0 / expected)
                        + ", " + nf1.format((oldmin - expected)));
                System.out.println(nf.format(d) + ", " + nf.format(min) + ", " + nf.format(expected) + ", "
                        + nf1.format(deviationMin * 100.0 / expected) + ", " + nf1.format(deviationMin));
            }
        }
    }

    /**
     * change in val at t by changing gseries coeff at midIdx
     */
    public static double[] changeToZeta(GSeries gSeries, final int initialPadding,
                                        double t, double initialValue, int midIdx, double increment) {
        double[] a0 = {0, 0};
        gSeries.incrementGValueAtIndex(midIdx, new double[]{increment, 0});
        
        double a00 = gSeries.evaluateZeta(t, initialPadding);
        a0[0] = (a00-initialValue)/ increment;
        
        gSeries.incrementGValueAtIndex(midIdx, new double[]{-increment, 0});
        
        // change second index

        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, increment});
        
        double a01 = gSeries.evaluateZeta(t, initialPadding);

        a0[1] = (a01-initialValue)/ increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, -increment});
        
        return a0;
    }
    
    public static double[] changeToDer(GSeries gSeries, final int initialPadding,
                                       double zero, double derAtZero, int midIdx, double increment) {
        double[] a0 = {0, 0};

        gSeries.incrementGValueAtIndex(midIdx, new double[]{increment, 0});
        double a0ValAtZero = gSeries.evalDer(
            zero, initialPadding, 0.00025* gSeries.spacing);
        a0[0] = (a0ValAtZero-derAtZero)/ increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{-increment, 0});
    
        // change second index
    
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, increment});
        a0ValAtZero = gSeries.evalDer(
              zero, initialPadding, 0.00025* gSeries.spacing);;
        a0[1] = (a0ValAtZero-derAtZero)/ increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, -increment});
    
        return a0;
    }
    
    public static double[] evalGSeriesIncrement(GSeries gSeries, int midIdx,
            final int initialPadding, double[] currentParams) {
        double zero = currentParams[0],  expectedDer = currentParams[1];
        double tmin = currentParams[3];
        double expectedMin = expectedDer<0?-currentParams[2]:currentParams[2];
        double evalAtMin = Interpolate.evaluateZeta(tmin, initialPadding, gSeries);

        double valAtZero = Interpolate.evaluateZeta(zero, initialPadding, gSeries);
        double derAtZero = Interpolate.evaluateDer(zero, initialPadding, gSeries);
        if(midIdx == -1){
           midIdx = gSeries.midIdx;
        }
        double[][] xy = new double[3][];
        xy[0] = changeToZeta(gSeries, initialPadding, zero, valAtZero, midIdx, 10);
        xy[1] = changeToDer(gSeries, initialPadding, zero, derAtZero, midIdx, 1);
        xy[2] = changeToZeta(gSeries, initialPadding, tmin, evalAtMin, midIdx, 10);
        
        double[] c = {-valAtZero, expectedDer-derAtZero, expectedMin-evalAtMin};
        if(Math.abs(c[0])<0.001 && Math.abs(c[1])<5){
            return new double[]{0, 0};
        }
        double sumxx = 0, sumxy = 0, sumyy = 0, sumcx = 0, sumcy = 0;
        for (int i = 0; i < xy.length; i++) {
            sumxx += xy[i][0]*xy[i][0];
            sumyy += xy[i][1]*xy[i][1];
            sumxy += xy[i][0]*xy[i][1];
            sumcx += c[i]*xy[i][0];
            sumcy += c[i]*xy[i][1];
        }
        double denom = sumxx*sumyy-sumxy*sumxy;
        // we will evaluate this as a simultaneous equation
        double[] g0incr = {
                (sumcx*sumyy-sumcy*sumxy)/(denom), 
                (sumcy*sumxx-sumcx*sumxy)/(denom)
                };
        if(!Double.isFinite(g0incr[0]) || !Double.isFinite(g0incr[1])){
            System.out.println(Arrays.toString(currentParams));
            System.out.println(Arrays.toString(c));
            System.out.println("valAtZero " + valAtZero + ", midIdx " + midIdx
                    + ", derAtZero " + derAtZero);
            System.out.println("denom " + denom);
            throw new IllegalStateException();
        }
        return g0incr;
    }

    public static double xmin(GSeries gSeries, final int initialPadding, double tmin) {
        double xa = tmin-0.01*gSeries.spacing;
        double xc = tmin+0.01*gSeries.spacing;
        double fa = Interpolate.evaluateZeta(xa, initialPadding, gSeries);
        double fb = Interpolate.evaluateZeta(tmin, initialPadding, gSeries);
        double fc = Interpolate.evaluateZeta(xc, initialPadding, gSeries);
        double r = (fa-fb)/(fb-fc);
        double xmin = (xa+tmin-r*(tmin+xc))/(2*(1-r));
        return xmin;
    }

}
