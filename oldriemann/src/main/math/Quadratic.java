/**
 * 
 */
package math;

import java.util.Arrays;

/**
 * coefficients of quadratic measured at equal increments
 * @author shankero
 *
 */
public class Quadratic {
    public static double[] coefficients(double t1, double incr, double[] z){
        double[] coefficients = new double[3];
        coefficients[2] = (z[0] + z[2] - 2*z[1])/(2*incr*incr);
        coefficients[1] = (z[2]-z[1])/incr - coefficients[2]*(2*t1+3*incr);
        coefficients[0] = ((t1+2*incr)*z[0]-t1*z[2])/(2*incr) + coefficients[2]*t1*(t1+2*incr);        
        return coefficients;
    }

    /**
     * @param args
     */
    public static void main(String[] args) {
        double[] z = new double[3];
        double incr = 1;
        double t1 = 10;
        for (int i = 0; i < z.length; i++) {
            double x  = t1 + i*incr;
            z[i] = 1.5 + x + 0.75*x*x;
            //z[i] = 500;
        }
        System.out.println("z " + Arrays.toString( z));
        double[] coeff = coefficients(t1, incr, z);
        System.out.println(Arrays.toString(coeff));
        for (int i = 0; i < z.length; i++) {
            double x  = t1+i*incr;
            double val = 0;
            double arg = 1;
            for (int j = 0; j < z.length; j++) {
                val += arg*coeff[j];
                arg *=x;
            }
            System.out.println(x + " " + val);
        }
    }

}
