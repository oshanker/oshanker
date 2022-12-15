package riemann;

import java.text.NumberFormat;

public class Poly7 {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }
    final double a, b, c;
    double d0, d1, d2;
    double t0, t1, t2;
    
    public Poly7(double a, double b, double c,
                 double a1, double b1, double c1) {
        this.a = a;
        this.b = b;
        this.c = c;
        d0 = a1;
        d1 = b1;
        d2 = c1;
        t0 = (a-b)*(a-c);
        t1 = (b-a)*(b-c);
        t2 = (c-a)*(c-b);
        t0 = d0/(t0*t0);
        t1 = d1/(t1*t1);
        t2 = d2/(t2*t2);
    }
    
    double der(double x) {
        if (Math.abs(x-a) < 1.0E-10 ) {
            return d0;
        }
        if (Math.abs(x-b) < 1.0E-10 ) {
            return d1;
        }
        if ( Math.abs(x-c) < 1.0E-10) {
            return d2;
        }
        double prod = (x-a)*(x-b)*(x-c);
        double ret = prod* (1/(x-a) + 1/(x-b) + 1/(x-c));
        double term = (t0*(x-b)*(x-c) + t1*(x-a)*(x-c)+ t2*(x-a)*(x-b));
        double dterm = (
            t0*(x-c) + t1*(x-c)+ t2*(x-b) +
            t0*(x-b) + t1*(x-a)+ t2*(x-a)
        );
        ret = ret*term + prod*dterm;
        return ret;
    }
    
    double eval(double x) {
        double prod = (x-a)*(x-b)*(x-c);
        double ret = prod*(t0*(x-b)*(x-c) + t1*(x-a)*(x-c)+ t2*(x-a)*(x-b));
        return ret;
    }
    
    double positionMax(double x0, double xa, double xb) {
        double der = der(x0);
        System.out.println(x0 + " " + der);
        if (Math.abs(der) < 1.0E-8) {
            return x0;
        }
        double incr = 0.001*(xb-xa);
        double next = x0 + incr;
        double derNext = der(next);
    
        for (int i = 0; i < 4; i++) {
            double oldx0 = x0;
            double oldder = der;
            der = derNext;
            x0 = next;
            next = oldx0 - oldder*(next-oldx0)/(derNext-oldder);
            derNext = der(next);
            System.out.println(next + " " + derNext);
            if (Math.abs(derNext) < 1.0E-8) {
                return next;
            }
        }
        return next;
    }
    
    public static void main(String[] args) {
        Poly7 poly7term = new Poly7(0, 1, 2, 2, -1, 2);
        tabulate(poly7term);
        double max1 = poly7term.positionMax(0.5, 0, 1);
        System.out.println("max1 " + nf.format(max1));
    }
    
    private static void tabulate(Poly7 poly7term) {
        for (double x = -0.1; x < 2.2; x += 0.05) {
            System.out.println(nf.format(x) +
                " " + nf.format(poly7term.eval(x)) +
                " " + nf.format(poly7term.der(x))
            );
        }
    }
}
