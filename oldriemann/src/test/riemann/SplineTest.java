package riemann;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Spline.Splinei;

public class SplineTest {

    @Test @Ignore
    public void testSplinei() {
        double xi = 1, h = 0.5;
        double yi = 2, yi1 = 2;
        double si = -1, si1 = 1;
        //2*(x-1)*(x-1.5) + 2
        Splinei splinei = new Splinei(xi, h, yi, yi1, si, si1);
        for (int i = 0; i < 3; i++) {
            double x = 1+0.25*i;
            double seconddDer = splinei.seconddDer(x);
            System.out.println(x + ", " + splinei.eval(x) + ", " + splinei.der(x)+ ", " + seconddDer);
            assertEquals(4.0, seconddDer, 0.0000001);
        }
        /**
         * h*(si2+4si1+si)=3*(yi2-yi)
         */
       }

        @Test
        public void testDerivativeMatching() {
            double xi = 0;
            double yi = 0, yi1 = 0;
            double slopei = 0, slopei1 = 64;
            //x^3-8x^2
            Splinei splinei = new Splinei(xi, 8, yi, yi1, slopei, slopei1);
            int N = 6;
            double[] si = new double[N];
            double[] y = new double[N];
            for (int i = 0; i < N; i++) {
                double x = 2*i;
                si[i] = splinei.der(x);
                y[i] = splinei.eval(x);
                System.out.println(i + ", " + y[i] + ", " + si[i] + ", " + splinei.seconddDer(x));
            }
            /**
             * h*(si2+4si1+si)=3*(yi2-yi)
             */
            double h = 2;
            double left = h*(si[3]+4*si[2]+si[1]);
            double right = 3*(y[3]-y[1]);
            assertEquals("From second der at first spline end :", left , right, 0.0000001 ) ;
            System.out.println("From second der at first spline end :" +  left + ", " + right ) ;
            left = h*(2*si[2]+4*si[1]);
            right = (5*y[2]-4*y[1]-y[0]);
            assertEquals("From first spline :", left , right, 0.0000001 ) ;
            System.out.println("From first spline :" +  left + ", " + right ) ;
            left = h*(2*si[N-3]+4*si[N-2]);
            right = (-5*y[N-3]+4*y[N-2]+y[N-1]);
            assertEquals("From last spline :", left , right, 0.0000001 ) ;
            System.out.println("From last spline :" +  left + ", " + right ) ;
        }

}
