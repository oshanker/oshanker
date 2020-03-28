package riemann;

import java.text.NumberFormat;
import java.util.Arrays;

public class NormalizedSpline extends BaseNormalizedSpline {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(6);
        nf.setMaximumFractionDigits(6);
        nf.setGroupingUsed(false);
    }

    /**
     * used only in test
     * @author shankero
     *
     */
    public static class Splinei {
        public final double xi, originalH;
        public final double yi, yi1;
        public final double si, si1;
        public Splinei(double xi, double h, double yi, double yi1,
                 double si, double si1) {
            this.xi = xi/h;
            this.yi = yi;
            this.yi1 = yi1;
            this.si = si;
            this.si1 = si1;
            this.originalH = h;
        }
        
        public double eval(double x){
            double ret = 0;
            double xifac = (x/originalH-xi);
            double xi1fac = (xifac-1);
            ret = xifac*xifac*(yi1 + (xi1fac)*(si1-2*yi1));
            ret += xi1fac*xi1fac*(yi + xifac*(si+2*yi));
            return ret;
        }
        
        public double der(double x){
            double ret = 0;
            double xifac = (x/originalH-xi);
            double xi1fac = (xifac-1);
            ret = 2*xifac*(yi1 + (xi1fac)*(si1-2*yi1));
            ret += xifac*xifac*((si1-2*yi1));
            ret += 2*xi1fac*(yi + (xifac)*(si+2*yi));
            ret += xi1fac*xi1fac*((si+2*yi));
            return ret;
        }

        public double seconddDer(double x){
            double ret = 0;
            double xifac = (x/originalH-xi);
            double xi1fac = (xifac-1);
            ret = 4*xifac*((si1-2*yi1));
            ret += 2*(yi1 + (xi1fac)*(si1-2*yi1));
            ret += 4*xi1fac*((si+2*yi));
            ret += 2*(yi + (xifac)*(si+2*yi));
            return ret;
        }
    }
    
    final double[] y;

    /**
     * x is from 0 to N-1
     * @param y
     */
    public NormalizedSpline(double[] y) {
        super(y.length);
        this.y = y;
        /**
         * si is being calculated
         */
       fit();
    }

    /**
     * populate diag and rhs with values needed by fitting procedure
     */
    protected void initSystem(double[] diag, double[] rhs) {
        diag[1] = 4.0;
        rhs[1] = (5*y[2]-4*y[1]-y[0]);
        for (int i = 2; i < diag.length-1; i++) {
            diag[i] = 4.0;
            rhs[i] = (3*(y[i+1]-y[i-1]));
        }
        diag[diag.length-1] = 4.0;
        rhs[diag.length-1] = (-5*y[N-3]+4*y[N-2]+y[N-1]);
    }

    /**
     * used in tests
     * @return
     */
    public double[] evalMid() {
        int index1 = 0;
        double[] rIm = new double[N-1];
        rIm[index1] = y[index1+2] - 3*(3*si[index1+1]+si[index1+2])/8;
        index1++;
        while(index1<N-2){
            rIm[index1] = (y[index1]+y[index1+1])/2 + (si[index1]-si[index1+1])/8;
            index1++;
        }
        rIm[index1] = y[index1-1] + 3*(3*si[index1]+si[index1-1])/8;
        return rIm;
    }

    public void evalMid(double[][] gSeries, int seriesOffset, int position) {
        int index1 = 0;
        gSeries[index1+seriesOffset][position] = y[index1+2] - 3*(3*si[index1+1]+si[index1+2])/8;
        index1++;
        while(index1<N-2){
            gSeries[index1+seriesOffset][position] = (y[index1]+y[index1+1])/2 + (si[index1]-si[index1+1])/8;
            index1++;
        }
        gSeries[index1+seriesOffset][position] = y[index1-1] + 3*(3*si[index1]+si[index1-1])/8;
    }


    public static void main(String[] args) {
        //x^3-8x^2
        int N = 7;
        double[] y = new double[N];
        double[] rIm = new double[N];
        double[] slope = new double[N];
        for (int i = 0; i < N; i++) {
            double x = 2*i;
            y[i] = f(x);
            slope[i] = 3*x*x - 16*x;
            rIm[i] = f(x+1);
        }
        System.out.println("in " + Arrays.toString(y));
        //rIm is f(x) at 1, 3, ..., 13
        System.out.println("expected " + Arrays.toString(rIm));
        // y is f(x) at 0, 2, ..., 12
        NormalizedSpline normalizedSpline = new NormalizedSpline(y);
        double[] actual = normalizedSpline.evalMid();
        System.out.println("actual " + Arrays.toString(actual));
        // slope is fprime(x) at 0, 2, ..., 12
        System.out.println("expected slopes " + Arrays.toString(slope));
        /**
         * si is double the slope.
         */
        System.out.println("actual slopes " + Arrays.toString(normalizedSpline.si));
    }

    private static double f(double x){
        return x*x*x - 8*x*x;
    }

}
