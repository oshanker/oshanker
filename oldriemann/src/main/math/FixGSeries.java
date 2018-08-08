package math;

import java.io.DataInputStream;
import java.io.File;
import java.text.NumberFormat;
import java.util.Arrays;

import riemann.Interpolate;

public class FixGSeries {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(8);
        nf.setMaximumFractionDigits(8);
        nf.setGroupingUsed(false);
    }

    public static void main(String[] args) throws Exception {
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

        int N = 1000102;
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

    public static double[] changeToZeta(GSeries gSeries, final int initialPadding, 
            double t, double valAtZero, int midIdx) {
        double[] a0 = {0,0};
        double increment = 10;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{increment, 0});
        double a0ValAtZero = Interpolate.evaluateZeta(t, initialPadding, gSeries);
        a0[0] = (a0ValAtZero-valAtZero)/increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{-increment, 0});
        
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, increment});
        a0ValAtZero = Interpolate.evaluateZeta(t, initialPadding, gSeries);
        a0[1] = (a0ValAtZero-valAtZero)/increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, -increment});
        return a0;
    }

    public static double[] changeToDer(GSeries gSeries, final int initialPadding, 
            double zero, double derAtZero, int midIdx) {
        double[] a0 = {0,0};
        double increment = 10;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{increment, 0});
        double a0ValAtZero = Interpolate.evaluateDer(zero, initialPadding, gSeries);
        a0[0] = (a0ValAtZero-derAtZero)/increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{-increment, 0});
        
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, increment});
        a0ValAtZero = Interpolate.evaluateDer(zero, initialPadding, gSeries);
        a0[1] = (a0ValAtZero-derAtZero)/increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, -increment});
        return a0;
    }

    public static double[] evalGSeriesIncrement(GSeries gSeries, int midIdx, 
            final int initialPadding, double[] currentParams) {
        double zero = currentParams[0],  expectedDer = currentParams[1];
        double tmin = currentParams[3];
        double expectedMin = expectedDer<0?-currentParams[2]:currentParams[2];
//        double xmin = FixGSeries.xmin(gSeries, initialPadding, tmin);
        double evalAtMin = Interpolate.evaluateZeta(tmin, initialPadding, gSeries);

        double valAtZero = Interpolate.evaluateZeta(zero, initialPadding, gSeries);
        double derAtZero = Interpolate.evaluateDer(zero, initialPadding, gSeries);
        if(midIdx == -1){
           midIdx = gSeries.midIdx;
        }
        double[][] xy = new double[3][];
        xy[0] = changeToZeta(gSeries, initialPadding, zero, valAtZero, midIdx);
        xy[1] = changeToDer(gSeries, initialPadding, zero, derAtZero, midIdx);
        xy[2] = changeToZeta(gSeries, initialPadding, tmin, evalAtMin, midIdx);
        
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
