package riemann;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Interpolate.Poly4;
import riemann.PolyInterpolate.BasePolySecondDer;
import riemann.PolyInterpolate.Poly5;

public class PolyInterpolateTest {

    @Test 
    public void testMultipleZeros() {
        //(x^2-3x+2)(1-(...))
        final double z0 = (3-Math.sqrt(5))/2, z1 = (3+Math.sqrt(5))/2;
        final double d0 = Math.sqrt(5), d1 = -Math.sqrt(5);
        double secondDer = -12;
        BasePolySecondDer basePolySecondDer = new BasePolySecondDer(
                z0, z1, d0, d1,  secondDer 
                ); 
        assertEquals(-5.0,basePolySecondDer.C, 0.000001);
        assertEquals(-12.0,basePolySecondDer.secondDerRHS(), 0.000001);
        assertEquals(-2.0, basePolySecondDer.eval(0), 0.000001);
        assertEquals(-1.0, basePolySecondDer.der(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.der(1.5), 0.000001);
        assertEquals(-0.3125,basePolySecondDer.eval(1.5), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval(2), 0.000001);
    }

    @Test
    public void testPoly5Interpolate() {
        double x0 = 9144.335987630979, x1 = 9144.449871243545;
        double d0 = 21.916679002964564, d1 = -18.66656999447713;
        double max = 0.8156186116467504;
        Poly5 poly5 = new Poly5(x0, x1, d0, d1, 74.30631140429409,
                max);
        System.out.println("second der at x1 " + poly5.secondDerRHS());
        System.out.println(poly5.positionMax + " " + poly5.eval(poly5.positionMax));
        System.out.println("poly5.C " + poly5.C);
        System.out.println("poly5.D " + poly5.D);
        Poly4 poly4 = new Poly4(x0, x1, d0, d1, 
                max);
        System.out.println("poly4 second der at x1 " + poly4.secondDerRHS());
        System.out.println(poly4.positionMax + " " + poly4.eval(poly4.positionMax));
        System.out.println("poly4.C " + poly4.C);
    }

    @Test
    public void testPoly5SecondDer() {
        double x0 = 2, x1 = 5;
        double xmin = Math.sqrt(1031.0/5);
        xmin = Math.sqrt(xmin);
        double minValue = Math.pow(xmin, 5) - 32 - (xmin-2)*1031;
        Poly5 poly5 = new Poly5(x0, x1, -951.0, 2094.0, 160,
                -minValue);
        assertEquals(2500.0, poly5.secondDerRHS(), 0.000001);
        assertEquals(xmin, poly5.positionMax, 0.000001);
        assertEquals(minValue, poly5.eval(poly5.positionMax), 0.000001);
        assertEquals(-1031, poly5.der(0), 0.000001);
        assertEquals(9.0, poly5.D, 0.000001);
    }
    
    @Test 
    public void testPoly5() {
        //(x^2-3x+2)(1-(...))
        Poly5 basePolySecondDer = new Poly5(1, 2, -1, 1, 0, 0.3125); 
        assertEquals(2.0, basePolySecondDer.eval1(0), 0.000001);
        assertEquals(-2.0, basePolySecondDer.eval(0), 0.000001);
        assertEquals(-1.0, basePolySecondDer.der(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.der(1.5), 0.000001);
        assertEquals(-0.3125,basePolySecondDer.eval(1.5), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval(1.5+Math.sqrt(1.25)), 0.000001);
        assertEquals(-1.0, basePolySecondDer.C, 0.000001);
        assertEquals(-0.0, basePolySecondDer.D, 0.000001);
    }

    @Test 
    public void testBasePolySecondDer() {
        //(x^2-3x+2)(1-(...))
        BasePolySecondDer basePolySecondDer = new BasePolySecondDer(1, 2, -1, 1, 0); 
        assertEquals(2.0, basePolySecondDer.eval1(0), 0.000001);
        assertEquals(-2.0, basePolySecondDer.eval(0), 0.000001);
        assertEquals(-1.0, basePolySecondDer.der(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.der(1.5), 0.000001);
        assertEquals(-0.3125,basePolySecondDer.eval(1.5), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval(1.5+Math.sqrt(1.25)), 0.000001);
        assertEquals(-1.0, basePolySecondDer.C, 0.000001);
    }

}
