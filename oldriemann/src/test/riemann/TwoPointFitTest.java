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
