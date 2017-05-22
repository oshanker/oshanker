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
        double[] begin = {//244.367502584863394599-delta, 244.367502584863394599+delta,
                244.920599505825861697-delta, 244.920599505825861697+delta};
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
        double[] roots = new double[]{243.8749480149, 244.158906912980683962,
                244.367502584863394599, 244.588579452072075626, 244.920599505825861697};
        double[] slopes = new double[]{17.618801272842216, -20.007038515089505,
                19.34333425934199, -27.175254718303503, 23.851008298572186};
        ZeroPoly zeroPoly = new ZeroPoly(roots, slopes);
        double t1 = 243.77756012466052947405878015472510;
        long offset = (long) 1.0E12;
        //double incr  = 2*Math.PI/(Math.log((offset+t1)/(2*Math.PI)));
        double incr = 0.24359904690399015;
        for (int i = 0; i < 5; i++) {
            double z = zeroPoly.eval(t1);
            double zetaFromRiemann = Riemann.riemann(t1, offset);
            System.out.println(t1 + " zetaFromRiemann " + zetaFromRiemann + " z "+ z);
            t1+=incr;
        }
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
