package riemann;

import static org.junit.Assert.*;

import java.text.NumberFormat;

import org.junit.Test;

import riemann.Interpolate.Poly3;
import riemann.Interpolate.Poly4;

public class InterpolateTest {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }

    @Test
    public void testMain() {
        final double z0 = 251.62429374748942, z1 = 252.0346709114582;
        final double d0 = 174.60414605716676, d1 = -207.33685149858513;
        final double max = 30.964994982487042;
        //244.2006260 zetaFromRiemann -0.7453399242908442
        //244.2627835 zetaFromRiemann -1.2321436376486554
        Poly3 poly = new Poly3(z0, z1, d0, d1);
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    );
        }
        System.out.println();
        System.out.println(nf.format(poly4.min) 
                + ", " + nf.format(poly4.der(poly4.min))
                + ", " + nf.format(poly4.eval(poly4.min))
        );
        double gram = 251.8163286404921;
        incr = (gram-251.62429374748942)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = 251.62429374748942 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    + ", poly " + nf.format(poly.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    );
        }
        assertEquals(30.597, poly4.eval(gram), 0.1);
    }

}
