package math;

import riemann.Interpolate;

public class FixGSeries {

    public static void main(String[] args) {
        // TODO Auto-generated method stub

    }

    public static double[] changeToZeta(GSeries gSeries, final int initialPadding, 
            double zero, double valAtZero, int midIdx) {
        double[] a0 = {0,0};
        double increment = 10;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{increment, 0});
        double a0ValAtZero = Interpolate.evaluateZeta(zero, initialPadding, gSeries);
        a0[0] = (a0ValAtZero-valAtZero)/increment;
        gSeries.incrementGValueAtIndex(midIdx, new double[]{-increment, 0});
        
        gSeries.incrementGValueAtIndex(midIdx, new double[]{0, increment});
        a0ValAtZero = Interpolate.evaluateZeta(zero, initialPadding, gSeries);
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
            final int initialPadding, double zero, double expectedDer) {
        double valAtZero = Interpolate.evaluateZeta(zero, initialPadding, gSeries);
        double derAtZero = Interpolate.evaluateDer(zero, initialPadding, gSeries);
        if(midIdx == -1){
           midIdx = gSeries.midIdx;
        }
        System.out.println(valAtZero + ", gSeries.midIdx " + midIdx
                + ", derAtZero " + derAtZero);
        double[][] xy = new double[2][];
        xy[0] = changeToZeta(gSeries, initialPadding, zero, valAtZero, midIdx);
        xy[1] = changeToDer(gSeries, initialPadding, zero, derAtZero, midIdx);
        double[] c = {-valAtZero, expectedDer-derAtZero};
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
