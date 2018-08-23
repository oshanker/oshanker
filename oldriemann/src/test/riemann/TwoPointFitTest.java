package riemann;

import static org.junit.Assert.*;

import org.junit.Test;

import riemann.Interpolate.Poly4;
import riemann.PolyInterpolate.BasePolySecondDer;
import riemann.PolyInterpolate.Poly5;

public class TwoPointFitTest {

    @Test
    public void testEval() {
        // x^3-8 - (x-2)*(125-8)/3 + x^4
        // 3x^2-39 + 4*x^3
        // 6x + 12*x^2
        TwoPointFit twoPointFit = new TwoPointFit(0, 70, -39, 0, 
                2, 0 + 16, -27 + 32, 12 + 48);
        assertEquals(625,twoPointFit.eval(5),0.000001);//625
        assertEquals(33,twoPointFit.eval(1),0.000001);//33
        assertEquals(0,twoPointFit.E + twoPointFit.F,0.000001);//33
    }

    @Test
    public void testPoly5SecondDer() {
        double x0 = 2, x1 = 5;
        double xmin = Math.sqrt(1031.0/5);
        xmin = Math.sqrt(xmin);
        double minValue = Math.pow(xmin, 5) - 32 - (xmin-2)*1031;
        final double d0 = -951.0;
        final double d1 = 2094.0;
        final int dd0 = 160;
        Poly5 poly5 = new Poly5(x0, x1, d0, d1, dd0,
                minValue);
        double x = -6.07828733661;
        System.out.println(x + ", (val) " + poly5.eval(x));
        final double positionMax = poly5.positionMax;
        assertEquals(xmin, positionMax, 0.000001);
        assertEquals(9.0, poly5.D, 0.000001);
        final double dd1 = 2500.0;
        assertEquals(dd1, poly5.secondDerRHS(), 0.000001);
        assertEquals(minValue, poly5.eval(positionMax), 0.000001);
        assertEquals(-1031, poly5.der(0), 0.000001);
        
        double mid = (x0+x1)/2;
        double midValue = poly5.eval(mid);
        TwoPointFit twoPointFit = new TwoPointFit(x0, 0, d0, dd0, x1, 0, d1, dd1);
        assertEquals(0d, twoPointFit.eval(-6.07828733661), 0.0000001);
        assertEquals(1,twoPointFit.E +  twoPointFit.F, 0.000001);
        assertEquals(minValue, twoPointFit.eval(positionMax), 0.000001);
        assertEquals(midValue, twoPointFit.eval(mid), 0.000001);
        assertEquals(midValue, TwoPointFit.mid(x1-x0, 0, d0, dd0,  0, d1, dd1), 
                0.000001);
    }

    @Test
    public void testInterpolate() {
        double z0 = -7, z1 = 2;
        double d0 = 108, d1 = -27;
        double max = 163.74433316206373;
        Poly4 poly = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin = poly.der(poly.positionMax);
        System.out.println(poly.C + ", " + poly.positionMax + ", " 
        + poly.eval(poly.positionMax)+ ", der " + derAtMin);
        assertEquals(163.74433316206373, poly.eval(poly.positionMax),0.000001);
        assertEquals(0d, Math.abs(derAtMin), 0.0000001);
        
        z0 = 2; z1 = 5;
        d0 = -27; d1 = 36;
        max = -23.74433316206373;
        Poly4 poly1 = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin1 = poly1.der(poly1.positionMax);
        System.out.println(poly1.C + ", " + poly1.positionMax + ", " 
        + poly1.eval(poly1.positionMax)+ ", der " + derAtMin1);
        assertEquals(-23.744333162063718, poly1.eval(poly1.positionMax),0.000001);
        assertEquals(0d, Math.abs(derAtMin1), 0.0000001);
        
        TwoPointFit twoPointFit = new TwoPointFit(
                0, poly.eval(-4.5), poly.der(-4.5), poly.secondDer(-4.5), 
                9, poly1.eval(4.5), poly1.der(4.5), poly1.secondDer(4.5));
        assertEquals(70, twoPointFit.eval(4.5),0.000001);
        double mid = TwoPointFit.mid(
                9, poly.eval(-4.5), poly.der(-4.5), poly.secondDer(-4.5), 
                 poly1.eval(4.5), poly1.der(4.5), poly1.secondDer(4.5));
        assertEquals(70, mid,0.000001);
    }
    
    @Test
    public void testEval1() {
        // x^3-8 - (x-2)*(125-8)/3 
        // 3x^2-39 
        // 6x 
        TwoPointFit twoPointFit = new TwoPointFit(0, 70, -39, 0, 
                2, 0, -27, 12);
        assertEquals(0,twoPointFit.eval(5),0.000001);//625
        assertEquals(0,twoPointFit.eval(-7),0.000001);//33
        assertEquals(1.0,twoPointFit.C + twoPointFit.D,0.000001);//33
        double xmax = Math.sqrt(13);
        assertEquals(163.74433316206373, twoPointFit.eval(-xmax),0.000001);
        assertEquals(-23.744333162063718, twoPointFit.eval(xmax),0.000001);
        System.out.println(-xmax + ", " + twoPointFit.eval(-xmax));
        System.out.println(xmax + ", " + twoPointFit.eval(xmax));
        
    }

}
