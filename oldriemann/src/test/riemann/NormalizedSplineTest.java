package riemann;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import riemann.NormalizedSpline.Splinei;

public class NormalizedSplineTest {
    private static double f(double x){
        return x*x*x - 8*x*x + 5*x - 40 ;
    }

    private static double der(double x){
        return 3*x*x - 16*x + 5;
    }


    @Test 
    public void testMid() {
        int N = 8;
        double[] x = new double[N];
        double[] y = new double[N];
        double[] midPointValues = new double[N];
        double[] slope = new double[N];
        final int mult = 3;
        for (int i = 0; i < N; i++) {
			x[i] = mult*i;
            y[i] = f(x[i]);
            midPointValues[i] = f(x[i]+mult/2.0);
            slope[i] = der(x[i]);
        }
        //x will be normalized from 0 to N-1
        NormalizedSpline normalizedSpline = new NormalizedSpline(y);
        assertEquals(slope.length-1, normalizedSpline.si.length);
        double[] actual = normalizedSpline.evalMid();
        for (int i = 0; i < N-1; i++) {
        	if(i>0) {
        		double testx =  mult/3.0;
                Spline.Splinei splinei = new Spline.Splinei (0, x[i+1]-x[i], y[i], y[i+1], 
                		slope[i], slope[i+1]);
                assertEquals("i = " + i + " :", f(x[i]+testx) , splinei.eval(testx), 0.0000001 ) ;
//                System.out.println("x "  + testx + " y " + f(x[i]+testx) + " pred " 
//                  + splinei.eval(testx));
                assertEquals("i = " + i + " :", slope[i] , normalizedSpline.si[i]/mult, 0.0000001 ) ;
        	}
            assertEquals(" i = " + i + " :", midPointValues[i] , actual[i], 0.0000001 ) ;
        }
        System.out.println("x " + Arrays.toString(x));
        System.out.println("y " + Arrays.toString(y));
        System.out.println("expected mid" + Arrays.toString(midPointValues));
        System.out.println("actual mid" + Arrays.toString(actual));
        System.out.println("expected slopes " + Arrays.toString(slope));
        System.out.println("actual slopes " + Arrays.toString(normalizedSpline.si));
     }
    
    @Test @Ignore
    public void test() {
    	//does this work?
        double xi = 0, x2 = 8;
        double yi = 0, y2 = 0;
        double slopei = 0, slope2 = 64*x2;
        //x^3-8x^2
        Splinei splinei = new Splinei(xi, x2, yi, y2, slopei, slope2);
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
               + ((si[i+1]+4*si[i]+si[i-1])/4) 
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
