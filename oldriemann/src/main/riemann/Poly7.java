package riemann;

import math.LinearEquation;

import java.text.NumberFormat;
import java.util.Arrays;

public class Poly7 implements Poly {
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
    double offset;
    
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
                 double a1, double b1, double c1,
                 double m0, double m1, double offset) {
        this(a, b, c, a1, b1, c1);
        this.m0 = m0;
        this.m1 = m1;
        this.offset = offset;
    
        setTermValues();
    }
    
    @Override
    public String toString() {
        return "Poly7{" +
            "a=" + a +
            ", b=" + b +
            ", c=" + c +
            ", d0=" + d0 +
            ", d1=" + d1 +
            ", d2=" + d2 +
            ", t0=" + t0 +
            ", t1=" + t1 +
            ", t2=" + t2 +
            ", m0=" + m0 +
            ", m1=" + m1 +
            ", offset=" + offset +
            '}';
    }
    
    public void setExtrema (double m0, double m1, double offset) {
        this.m0 = m0;
        this.m1 = m1;
        this.offset = offset;
    }
    
    public double setTermValues() {
        double positionMax0 = positionMax((a + b) / 2, a, b);
        double currentMax0 = eval(positionMax0);
        double positionMax1 = positionMax((b + c) / 2, b, c);
        double currentMax1 = eval(positionMax1);
        int iter = 0;
        double deviation = (Math.abs(m0 -currentMax0) + Math.abs(m1 -currentMax1));
        System.out.println("input: deviation " + deviation);
        while(deviation>1.0E-8) {
            double[][] coeff = new double[2][2];
            populateCoeff(currentMax0, currentMax1, 0, coeff);
            populateCoeff(currentMax0, currentMax1, 1, coeff);
            LinearEquation linearEquation = new LinearEquation(coeff);
            double[] neededZetaIncrement = {
                (m0 -currentMax0),
                (m1 -currentMax1)
            };
            double[] solution = linearEquation.solve(
                neededZetaIncrement
            );
            System.out.println(iter++ + " Required increment " );
            System.out.println( Arrays.toString(solution));
            incrementTermTemp(solution);
            positionMax0 = positionMax(positionMax0, a, b);
            currentMax0 = eval(positionMax0);
            positionMax1 = positionMax(positionMax1, b, c);
            currentMax1 = eval(positionMax1);
            deviation = (Math.abs(m0 -currentMax0) + Math.abs(m1 -currentMax1));
            System.out.println("deviation " + deviation);
            if (iter>8) {
                break;
            }
        }
        return deviation;
    }
    
    private void populateCoeff( double currentMax0, double currentMax1, int idx, double[][] coeff) {
    
        double incr = 0.1;
        double[] a1b1 = {0, 0};
        a1b1[idx] = incr;
        incrementTermTemp(a1b1);
        double pNextMax0 = positionMax((a + b) / 2, a, b);
        double cNextMax0 = eval(pNextMax0);
        double pNextMax1 = positionMax((b + c) / 2, b, c);
        double cNextMax1 = eval(pNextMax1);
        coeff[0][idx] = (cNextMax0- currentMax0)/ incr;
        coeff[1][idx] = (cNextMax1- currentMax1)/ incr;
        unsetTerm();
    }
    
    void setTerm(double a1, double b1, double offset) {
        poly7term = new Poly7term(a, b, c, a1, b1, offset);
    }
    
    void incrementTermTemp(double[] a1b1) {
        oldTerm = poly7term;
        if(poly7term==null) {
            poly7term = new Poly7term(a, b, c, a1b1[0], a1b1[1], offset);
        } else {
            poly7term = new Poly7term(a, b, c,
                poly7term.A + a1b1[0],
                poly7term.B + a1b1[1], offset);
        }
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
    
    @Override
    public double der(double x) {
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
    
    @Override
    public double eval(double x) {
        double prod = (x-a)*(x-b)*(x-c);
        double ret = prod*(t0*(x-b)*(x-c) + t1*(x-a)*(x-c)+ t2*(x-a)*(x-b));
        if(poly7term != null) {
            ret += poly7term.eval(x);
        }
        return ret;
    }
    
    @Override
    public double getPositionMax() {
        return positionMax((a + b) / 2, a, b);
    }
    
    @Override
    public double secondDer(double t) {
        double incr = 0.001*(b-a);
        double ret = (der(t + incr) - der(t - incr)) / (2 * incr);
        return ret;
    }
    
    double positionMax(double x0, double xa, double xb) {
        double derx0 = 0;
        while (xb-xa>0.001) {
            double signumxa = Math.signum(der(xa));
            if (signumxa == 0) {
                return xa;
            }
            double signumxb = Math.signum(der(xb));
            if (signumxb == 0) {
                return xb;
            }
            derx0 = der(x0);
            double signumx0 = Math.signum(derx0);
            if (signumx0 == 0) {
                return x0;
            }
            if(signumx0 == signumxa) {
                xa = x0;
            } else {
                xb = x0;
            }
            x0 = (xa+xb)/2;
        }
//        System.out.println("x0 " + x0 + " xa " + xa + " xb " + xb
//        + " derx0 " + derx0);
        
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
            if (next < xa || next > xb){
                throw new IllegalStateException(xa + " next " + next
                + " xb " + xb);
            }
            derNext = der(next);
//            System.out.println(next + " derNext " + derNext);
            if (Math.abs(derNext) < 1.0E-8) {
                return next;
            }
        }
        if(!Double.isFinite(next)){
            throw new IllegalStateException(xa + " ");
        }
        return next;
    }
    
    public static void main(String[] args) {
        //B = 0.5
        Poly7 poly7 = new Poly7(0, 1, 2, 2, -1, 2,
            0.4589742535338246, -0.31082610538567645, 0);
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        System.out.println("============= " );
        //A = 0.5
        poly7 = new Poly7(0, 1, 2, 2, -1, 2,
            0.41689197105413617, -0.2700198241673988, 0);
        System.out.println("max0 " + poly7.evalMax0());
        System.out.println("max1 " + poly7.evalMax1());
        //tabulate(poly7);
        //y = y0 + (y1-y0)*(x-x0)/(x1-x0)
        //y0*(x1-x0) + (y1-y0)*(x-x0)
        //x = x0-y0*(x1-x0)/(y1-y0)
    }
    
    public void tabulate(double xa, double xb, int steps) {
        double incr = (xb-xa)/(steps-1);
        
        for (int i = 0; i < steps; i++) {
            double x = xa + i*incr;
            System.out.println(nf.format(x) +
                " " + nf.format(eval(x)) +
                " der " + nf.format(der(x))
            );
        }
        
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
