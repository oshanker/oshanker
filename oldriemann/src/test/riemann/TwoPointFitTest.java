package riemann;

import static org.junit.Assert.*;

import org.junit.Test;

import riemann.Interpolate.Poly4;

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
    public void testInterpolate() {
        double z0 = -7, z1 = 2;
        double d0 = 108, d1 = -27;
        double max = 163.74433316206373;
        Poly4 poly = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin = poly.der(poly.positionMax);
        System.out.println(poly.C + ", " + poly.positionMax + ", " 
        + poly.eval(poly.positionMax)+ ", der " + derAtMin);
        assertEquals(163.74433316206373, poly.eval(poly.positionMax),0.000001);
        assertEquals(0d, derAtMin, 0.0000001);
        
        z0 = 2; z1 = 5;
        d0 = -27; d1 = 36;
        max = -23.74433316206373;
        Poly4 poly1 = new Poly4(z0, z1, d0, d1, max);
        final double derAtMin1 = poly1.der(poly1.positionMax);
        System.out.println(poly1.C + ", " + poly1.positionMax + ", " 
        + poly1.eval(poly1.positionMax)+ ", der " + derAtMin1);
        assertEquals(-23.744333162063718, poly1.eval(poly1.positionMax),0.000001);
        assertEquals(0d, derAtMin1, 0.0000001);
        
        TwoPointFit twoPointFit = new TwoPointFit(
                0, poly.eval(-4.5), poly.der(-4.5), poly.secondDer(-4.5), 
                9, poly1.eval(4.5), poly1.der(4.5), poly1.secondDer(4.5));
        assertEquals(70, twoPointFit.eval(4.5),0.000001);
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
