package riemann;

import java.text.NumberFormat;

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
    };

    public static void main(String[] args) {
        final double z0 = 244.158906912980683962, z1 = 244.367502584863394599;
        final double d0 = -20.007604626096071598, d1 = 19.343950349024609636;
        //244.2006260 zetaFromRiemann -0.7453399242908442
        //244.2627835 zetaFromRiemann -1.2321436376486554
        Poly3 poly = new Poly3(z0, z1, d0, d1);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) + ", " + nf.format(poly.eval(x)));
        }
        double xmin = z0 - d0*(z1-z0)/(d1-d0);
        System.out.println(nf.format(xmin) + ", " + nf.format(poly.eval(xmin)));

    }

}
