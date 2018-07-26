package riemann;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import riemann.PolyInterpolate.BasePolySecondDer;
import riemann.PolyInterpolate.Poly5;

public class PolyInterpolateTest {

    @Test
    public void testPoly5() {
        //(x^2-3x+2)(1-(...))
        Poly5 basePolySecondDer = new Poly5(1, 2, -1, 1, 0, 0.3125); 
        assertEquals(2.0, basePolySecondDer.eval1(0), 0.000001);
        assertEquals(-2.0, basePolySecondDer.eval2(0), 0.000001);
        assertEquals(-1.0, basePolySecondDer.der(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.der(1.5), 0.000001);
        assertEquals(-0.3125,basePolySecondDer.eval2(1.5), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval2(1.5+Math.sqrt(1.25)), 0.000001);
        assertEquals(-1.0, basePolySecondDer.C, 0.000001);
        assertEquals(-0.0, basePolySecondDer.D, 0.000001);
    }

    @Test @Ignore
    public void testBasePolySecondDer() {
        //(x^2-3x+2)(1-(...))
        BasePolySecondDer basePolySecondDer = new BasePolySecondDer(1, 2, -1, 1, 0); 
        assertEquals(2.0, basePolySecondDer.eval1(0), 0.000001);
        assertEquals(-2.0, basePolySecondDer.eval2(0), 0.000001);
        assertEquals(-1.0, basePolySecondDer.der(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.der(1.5), 0.000001);
        assertEquals(-0.3125,basePolySecondDer.eval2(1.5), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval2(1.5+Math.sqrt(1.25)), 0.000001);
        System.out.println(basePolySecondDer.C);
    }

}
