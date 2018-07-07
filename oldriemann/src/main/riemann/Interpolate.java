package riemann;

import java.text.NumberFormat;
import java.util.Arrays;

public class Interpolate {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }
    
    public static class Poly3{
        final double z0, z1;
        final double d0, d1;
        double denom;
        public Poly3(double z0, double z1, double d0, double d1) {
            this.z0 = z0;
            this.z1 = z1;
            this.d0 = d0;
            this.d1 = d1;
            denom = (z0-z1)*(z0-z1);
        }
        public double eval(double x){
            double ret = (x-z0)*(x-z1)/denom;
            ret *= ((x-z0)*d1 + (x-z1)*d0);
            return ret;
        }
        public double der(double x){
            double ret = (
                    (d1+d0)*(x-z0)*(x-z1) + (x-z0)*((x-z0)*d1 + (x-z1)*d0) + (x-z1)*((x-z0)*d1 + (x-z1)*d0)
                    )/denom;
            return ret;
            
        }
    };

    public static class Poly4 extends Poly3{
        double C;
        double min;
        public Poly4(double z0, double z1, double d0, double d1, double max) {
            super(z0, z1, d0, d1);
            if(d0<0){max = -max;}
            double[] oldest = new double[]{z0, d0, 0};
            double[] upper = new double[]{z1, d1, 1};
            double[] wts = new double[]{Math.abs(d1),Math.abs(d0)};
            double xmin = z0 - d0*(z1-z0)/(d1-d0);
            //C = estimateC( max, est[0]);
            for (int i = 0; i < 4; i++) {
               C = estimateC( max, xmin);
               xmin(oldest, upper, wts);
               System.out.println(Arrays.toString(oldest)+ ", " + Arrays.toString(upper));
               if(oldest[2]> 99){
                   xmin = oldest[0];
                   break;
               }
               wts[0] = Math.abs(upper[1]);
               wts[1] = Math.abs(oldest[1]);
               if(wts[0]>wts[1]){
                   xmin = oldest[0];
               }else {
                   xmin = upper[0];
               }
            }
            C = estimateC( max, xmin);
            min = xmin;
            System.out.println(C + ", " + xmin + ", " + eval(xmin)+ ", der " + der(xmin));
        }

        private double estimateC( double max, double xmin) {
            double mult = (xmin-z0)*(xmin-z1);
            mult = mult*mult/denom;
            return (max - super.eval(xmin))/mult;
        } 
        
        void xmin(double[] oldest, double[] upper, double[] wts){
            double x0 = oldest[0];
            double x1 = upper[0];
            double xmin = (wts[0]*x0 + wts[1]*x1)/(wts[0]+wts[1]);
            double dermin = der(xmin);
            if(Math.abs(dermin)<0.01){
                oldest[0]=xmin;
                oldest[1]=dermin;
                oldest[2]=100;
            } else if(Math.signum(dermin) == Math.signum(oldest[1])){
                oldest[0]=xmin;
                oldest[1]=dermin;
            } else {
                upper[0]=xmin;
                upper[1]=dermin;
            }
        }
        public double eval(double x){
            double mult = (x-z0)*(x-z1);
            mult = mult*mult/denom;
            double ret = super.eval(x)+C*mult;
            return ret;
        }        
        public double der(double x){
            double ret = super.der(x);
            ret += 2*C*(x-z0)*(x-z1)*(2*x-z1-z0)/denom;
            return ret;
        }
    }
    
    public static void main(String[] args) {
        testInitialZeros();
        //debug();

    }

    private static void debug() {
        final double z0 = 1, z1 = 2;
        final double d0 = -1, d1 = 1;
        final double max = 0.25;
        Poly3 poly = new Poly3(z0, z1, d0, d1);
        int N = 11;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    );
        }
    }

    private static void testInitialZeros() {
        final double z0 = 244.158906912980683962, z1 = 244.367502584863394599;
        final double d0 = -20.007604626096071598, d1 = 19.343950349024609636;
        final double max = 1.232146174810101691;
        //244.2006260 zetaFromRiemann -0.7453399242908442
        //244.2627835 zetaFromRiemann -1.2321436376486554
        Poly3 poly = new Poly3(z0, z1, d0, d1);
        double xmin = 244.2627835;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly.eval(x))
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    + ", " + nf.format(poly4.der(x))
                    );
        }
        System.out.println(nf.format(xmin) 
                + ", " + nf.format(poly.der(xmin))
                + ", " + nf.format(poly4.der(xmin))
                );
        System.out.println(nf.format(xmin) 
                + ", " + nf.format(poly.eval(xmin))
                + ", " + nf.format(poly4.eval(xmin))
        );
        System.out.println(nf.format(poly4.min) 
                + ", " + nf.format(poly4.der(poly4.min))
                + ", " + nf.format(poly4.eval(poly4.min))
        );
        incr = (244.2627835-244.26257831011117)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = 244.26257831011117 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    + ", poly " + nf.format(poly.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    );
        }
    }

}
