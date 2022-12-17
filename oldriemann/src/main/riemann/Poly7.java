package riemann;

import java.text.NumberFormat;

public class Poly7 {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }
    Poly7term poly7term = null;
    final double a, b, c;
    final double d0, d1, d2;
    double t0, t1, t2;
    double m0, m1;
    private Poly7term oldTerm;
    
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
    
    public Poly7(double a, double b, double c,
                 double a1, double b1, double c1, double m0, double m1) {
        this(a, b, c, a1, b1, c1);
        this.m0 = m0;
        this.m1 = m1;
        double positionMax0 = positionMax((a + b) / 2, a, b);
        double currentMax0 = eval(positionMax0);
        double positionMax1 = positionMax((b + c) / 2, b, c);
        double currentMax1 = eval(positionMax1);
        if(Math.abs(m0-currentMax0) + Math.abs(m1-currentMax1)>1.0E-8) {
            double[][] coeff = new double[2][2];
            populateCoeff(currentMax0, currentMax1, 1, coeff);
            double delB0 = (m0-currentMax0)/coeff[0][1];
            System.out.println(delB0);
            double delB1 = (m1-currentMax1)/coeff[0][1];
            System.out.println(delB1);
    
            poly7term = new Poly7term(a, b, c,0, (delB0+delB1)/2);
        }
        
    }
    
    private double populateCoeff( double currentMax0, double currentMax1, int idx, double[][] coeff) {
        double deviation = (Math.abs(m0-currentMax0) + Math.abs(m1-currentMax1));
        if(deviation<1.0E-8){
            return deviation;
        }
        double incr = 0.1;
        double[] a1b1 = {0, 0};
        a1b1[idx] = incr;
        setTermTemp(a1b1);
        double pNextMax0 = positionMax((a + b) / 2, a, b);
        double cNextMax0 = eval(pNextMax0);
        double pNextMax1 = positionMax((b + c) / 2, b, c);
        double cNextMax1 = eval(pNextMax1);
        coeff[0][idx] = (cNextMax0- currentMax0)/ incr;
        coeff[1][idx] = (cNextMax1- currentMax1)/ incr;
        unsetTerm();
        return deviation;
    }
    
    void setTerm(double a1, double b1) {
        poly7term = new Poly7term(a, b, c, a1, b1);
    }
    
    void setTermTemp(double[] a1b1) {
        oldTerm = poly7term;
        poly7term = new Poly7term(a, b, c, a1b1[0], a1b1[1]);
    }
    
    void unsetTerm() {
        poly7term = oldTerm;
        oldTerm = null;
    }
    
    double evalMax0() {
        double pmax0 = positionMax((a+b)/2, a, b);
        double max0 = eval(pmax0);
        if(Math.abs(max0-m0) < 1.0E-8) {
            System.out.println("m0 OK");
        } else {
            System.out.println("m0 dev " + (max0-m0));
        }
        return max0;
    }
    
    double evalMax1() {
        double pmax1 = positionMax((c+b)/2, b, c);
        double max1 = eval(pmax1);
        if(Math.abs(max1-m1) < 1.0E-8) {
            System.out.println("m1 OK");
        } else {
            System.out.println("m1 dev " + (max1-m1));
        }
        return max1;
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
        if(poly7term != null) {
            ret += poly7term.der(x);
        }
        return ret;
    }
    
    double eval(double x) {
        double prod = (x-a)*(x-b)*(x-c);
        double ret = prod*(t0*(x-b)*(x-c) + t1*(x-a)*(x-c)+ t2*(x-a)*(x-b));
        if(poly7term != null) {
            ret += poly7term.eval(x);
        }
        return ret;
    }
    
    double positionMax(double x0, double xa, double xb) {
        double der = der(x0);
        if (Math.abs(der) < 1.0E-8) {
            return x0;
        }
        double incr = 0.001*(xb-xa);
        double next = x0 + incr;
        double derNext = der(next);
        if (Math.abs(derNext) < 1.0E-8) {
            return next;
        }
    
        for (int i = 0; i < 15; i++) {
            double oldx0 = x0;
            double oldder = der;
            der = derNext;
            x0 = next;
            next = oldx0 - oldder*(next-oldx0)/(derNext-oldder);
            derNext = der(next);
            if (Math.abs(derNext) < 1.0E-8) {
                return next;
            }
        }
        return next;
    }
    
    public static void main(String[] args) {
        //B = 0.5
        Poly7 poly7 = new Poly7(0, 1, 2, 2, -1, 2,
            0.4589742535338246, -0.31082610538567645);
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        //A = 0.5
        poly7 = new Poly7(0, 1, 2, 2, -1, 2,
            0.41689197105413617, -0.2700198241673988);
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        //tabulate(poly7);
        //y = y0 + (y1-y0)*(x-x0)/(x1-x0)
        //y0*(x1-x0) + (y1-y0)*(x-x0)
        //x = x0-y0*(x1-x0)/(y1-y0)
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
