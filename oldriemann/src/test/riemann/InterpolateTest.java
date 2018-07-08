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
        final double z0 = 244.158906912980683962, z1 = 244.367502584863394599;
        final double d0 = -20.007604626096071598, d1 = 19.343950349024609636;
        final double max = 1.232146174810101691;
        //244.2006260 zetaFromRiemann -0.7453399242908442
        //244.2627835 zetaFromRiemann -1.2321436376486554
        Poly3 poly = new Poly3(z0, z1, d0, d1);
        double xmin = 244.2627835;
        Poly4 poly4 = new Poly4(z0, z1, d0, d1, max);
        int N = 6;
        double incr = (z1-z0)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = z0 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly.eval(x))
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    + ", " + nf.format(poly4.der(x))
                    );
        }
        System.out.println(nf.format(xmin) 
                + ", " + nf.format(poly.der(xmin))
                + ", " + nf.format(poly4.der(xmin))
                );
        System.out.println(nf.format(xmin) 
                + ", " + nf.format(poly.eval(xmin))
                + ", " + nf.format(poly4.eval(xmin))
        );
        System.out.println(nf.format(poly4.min) 
                + ", " + nf.format(poly4.der(poly4.min))
                + ", " + nf.format(poly4.eval(poly4.min))
        );
        incr = (244.2627835-244.26257831011117)/(N-1);
        for (int i = 0; i < N; i++) {
            double x = 244.26257831011117 + i*incr;
            System.out.println(nf.format(x) 
                    + ", " + nf.format(poly4.eval(x))
                    + ", der " + nf.format(poly4.der(x))
                    + ", poly " + nf.format(poly.eval(x))
                    + ", der " + nf.format(poly.der(x))
                    );
        }
    }

}
