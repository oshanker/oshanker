package riemann;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Spline.Splinei;

public class SplineTest {

    @Test //@Ignore
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
    
    @Test //@Ignore
    public void test1() {
        double xi = 0, x2 = 3;
        double yi = -40, y2 = -70;
        double slopei = 5, slope2 = -16;
        //47.125
        Splinei splinei = new Splinei(xi, x2, yi, y2, slopei, slope2);
        assertEquals(-47.125, splinei.eval(1.5), 0.00001);
    }   
    
    @Test //@Ignore
    public void test2() {
        double xi = 3, x2 = 6;
        double yi = -70, y2 = -82;
        double slopei = -16, slope2 = 17;
        //47.125
        Splinei splinei = new Splinei(xi, x2-xi, yi, y2, slopei, slope2);
        assertEquals(-88.375, splinei.eval(xi+1.5), 0.00001);
    }   
	    @Test //@Ignore
	    public void testSplinei2() {
	        double xi = 0, xmax = 10;
	        double yi = 0, yi1 = 200;
	        double slopei = 0, slopei1 = 140;
	        //x^3-8x^2
	        //3*x*x-16*x
	        Splinei splinei = new Splinei(xi, xmax, yi, yi1, slopei, slopei1);
	        int N = 5;
	        double[] si = new double[N];
	        double[] y = new double[N];
	        for (int i = 0; i < N; i++) {
	            double x = 2*i;
	            si[i] = splinei.der(x);
	            assertEquals(3*x*x-16*x, si[i], 0.0001);
	            y[i] = splinei.eval(x);
	            System.out.println(i + ", " + x + ", " + y[i] + ", " 
	            + si[i] + ", " + splinei.seconddDer(x));
	        }
	        /**
	         * h*(si2+4si1+si)=3*(yi2-yi)
	         */
	        System.out.println("From second der at x = 4 :" + (2*(si[3]+4*si[2]+si[1])) + ", " + (3*(y[3]-y[1])));
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
