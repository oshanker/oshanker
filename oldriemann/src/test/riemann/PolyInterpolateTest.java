package riemann;

import static org.junit.Assert.*;

import org.junit.Test;

import riemann.PolyInterpolate.BasePolySecondDer;

public class PolyInterpolateTest {

    @Test
    public void testMain() {
        BasePolySecondDer basePolySecondDer = new BasePolySecondDer(1, 2, -1, 1, 0); 
        assertEquals(2.0, basePolySecondDer.eval1(0), 0.000001);
        assertEquals(-2.0, basePolySecondDer.eval2(0), 0.000001);
        assertEquals(-1.0, basePolySecondDer.der(1), 0.000001);
        assertEquals(0.0, basePolySecondDer.eval2(1.5+Math.sqrt(1.25)), 0.000001);
    }

}
