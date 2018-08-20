package riemann;

import static org.junit.Assert.*;

import java.text.NumberFormat;
import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Interpolate.Poly3;
import riemann.Interpolate.Poly4;

public class InterpolateTest {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }
    
    public double eval(double x, Poly3 poly, double C){
        double mult = (x-poly.z0)*(x-poly.z1);
        mult = mult*mult/poly.denom;
        double ret = poly.eval1(x)+C*mult;
        return ret;
    }        

    @Test @Ignore
    public void testSlowGrowth() {
        final double z0 = 247.1270084390591, z1 = 247.25934052902014;
        final double d0 = -5.46251113952012, d1 = 4.963052066753181;
        final double max = 0.19014781349558818;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        double gram = 247.18794676831632;
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        System.out.println(nf.format(gram) 
                + ", " + nf.format(poly4.eval(gram))
                + ", der " + nf.format(poly4.der(gram))
                );
        assertEquals(-0.189428825055622, poly4.eval(gram), 0.00005);
        System.out.println();
    }
    
    @Test 
    public  void debug() {
        final double z0 = 1, z1 = 2;
        final double d0 = -1, d1 = 1;
        final double max = -0.25;
        Poly4 poly = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin = poly.der(poly.positionMax);
        System.out.println(poly.C + ", " + poly.positionMax + ", " 
        + poly.eval(poly.positionMax)+ ", der " + derAtMin);
        assertEquals(2.0d, poly.secondDerRHS(), 0.0000001);
        assertEquals(2.0d, poly.secondDer(1.75), 0.0000001);
        assertEquals(0d, poly.C, 0.0000001);
        assertEquals(0d, derAtMin, 0.0000001);
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

    @Test @Ignore
    public void testConvergence() {
        final double z0 = 251.306919355202, z1 = 251.62429374748942;
        final double d0 = -51.032529476335846, d1 = 174.60414605716676;
        final double max = 11.500212986748618;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        double gram = 251.3291305486841;
        incr = (gram-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        assertEquals(-1.283921548, poly4.eval(gram), 0.1);
        System.out.println();
        gram = 251.57272959458808;
        incr = (gram-poly4.positionMax)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = poly4.positionMax + i*incr;
            System.out.println("n 33, " + nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        assertEquals(-7.45540885886704, poly4.eval(gram), 0.2);
    }

    @Test @Ignore
    public void testMain() {
        final double z0 = 251.62429374748942, z1 = 252.0346709114582;
        final double d0 = 174.60414605716676, d1 = -207.33685149858513;
        final double max = 30.964994982487042;
        //244.2006260 zetaFromRiemann -0.7453399242908442
        //244.2627835 zetaFromRiemann -1.2321436376486554
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        System.out.println();
        System.out.println(nf.format(poly4.positionMax) 
                + ", " + nf.format(poly4.der(poly4.positionMax))
                + ", " + nf.format(poly4.eval(poly4.positionMax))
        );
        double gram = 251.8163286404921;
        incr = (gram-251.62429374748942)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = 251.62429374748942 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        assertEquals(30.597, poly4.eval(gram), 0.1);
    }

}
