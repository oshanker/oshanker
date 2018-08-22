package riemann;

import static org.junit.Assert.*;

import org.junit.Test;

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

}
