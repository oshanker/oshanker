package math;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Test;

public class ZeroPolyTest {

    @Test
    public void test() {
        double[] roots = new double[]{-1, 1};
        double[] slopes = new double[]{-4, 4};
        ZeroPoly zeroPoly = new ZeroPoly(roots, slopes);
        double t1 = -0.5;
        double incr = 1;
        double[] coeff = zeroPoly.coefficients(t1 , incr );
        //System.out.println("coeff " + Arrays.toString(coeff));
        double[] expected = new double[]{-2, 0,2};
        assertArrayEquals(expected, coeff, 1.0E-10);
    }

}
