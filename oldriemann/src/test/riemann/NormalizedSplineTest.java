package riemann;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import riemann.NormalizedSpline.Splinei;

public class NormalizedSplineTest {
    private static double f(double x){
        return x*x*x - 8*x*x + 5*x - 40 ;
    }


    @Test
    public void testMid() {
        int N = 8;
        double[] y = new double[N];
        double[] rIm = new double[N];
        for (int i = 0; i < N; i++) {
            double x = 2*i;
            y[i] = f(x);
            rIm[i] = f(x+1);
        }
        NormalizedSpline normalizedSpline = new NormalizedSpline(y);
        double[] actual = normalizedSpline.evalMid();
        for (int i = 0; i < N-1; i++) {
            assertEquals("predicted i = " + i + " :", rIm[i] , actual[i], 0.0000001 ) ;
        }
    }
    
    @Test @Ignore
    public void test() {
        double xi = 0, h = 8;
        double yi = 0, yi1 = 0;
        double slopei = 0, slopei1 = 64*h;
        //x^3-8x^2
        Splinei splinei = new Splinei(xi, h, yi, yi1, slopei, slopei1);
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
        for (int i = 2; i < N-1; i++) {
            System.out.println("From second der at i = " + i + " :" 
               + ((si[i+1]+4*si[i]+si[i-1])/(N-1)) 
                    + ", " + (3*(y[i+1]-y[i-1])));
            assertEquals("From second der at i = " + i 
                    + " :", ((si[i+1]+4*si[i]+si[i-1])/4) 
                    ,  (3*(y[i+1]-y[i-1])), 0.000001);
        }
        double left = (2*si[2]+4*si[1])/4;
        double right = (5*y[2]-4*y[1]-y[0]);
        System.out.println("From first spline :" +  left + ", " + right ) ;
        assertEquals("From first spline :", left , right, 0.0000001 ) ;
        left = (2*si[N-3]+4*si[N-2])/4;
        right = (-5*y[N-3]+4*y[N-2]+y[N-1]);
        System.out.println("From last spline :" +  left + ", " + right ) ;
        assertEquals("From last spline :", left , right, 0.0000001 ) ;
   }

}
