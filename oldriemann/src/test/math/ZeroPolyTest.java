package math;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Riemann;

public class ZeroPolyTest {
    @Test @Ignore
    public void riemann() {
        double delta = 0.001;
        double[] begin = {243.8749480149-delta, 243.8749480149+delta,
                244.158906912980683962-delta, 244.158906912980683962+delta};
        double[] zetaFromRiemann = new double[begin.length];
        long offset = (long) 1.0E12;
        for (int i = 0; i < begin.length; i++) {
            zetaFromRiemann[i] = Riemann.riemann(begin[i], offset);
            double der = 0;
            if(i%2==1){
                der = (zetaFromRiemann[i]-zetaFromRiemann[i-1])/(2*delta);
                System.out.println("der " + der);
            }
            System.out.println("zetaFromRiemann " + zetaFromRiemann[i]);
        }
    }
    
    @Test 
    public void root() {
        double[] roots = new double[]{243.8749480149, 244.158906912980683962};
        double[] slopes = new double[]{17.618801272842216, -20.007038515089505};
        ZeroPoly zeroPoly = new ZeroPoly(roots, slopes);
        double t1 = 243.77756012466052947405878015472510;
        double z = zeroPoly.eval(t1);
        System.out.println("z "+ z);
    }
    
    @Test @Ignore
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
