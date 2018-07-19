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
    final double[] coefficients;
    public Quadratic(double firstPoint, double incr, double[] values){
       this.coefficients = coefficients(firstPoint, incr, values);
    }
    
    public double eval(double x ){
        double val = 0;
        double arg = 1;
        for (int j = 0; j < 3; j++) {
            val += arg*coefficients[j];
            arg *=x;
        }
        return val;
    }
 
    public double der(double x ){
        return coefficients[1] + 2* coefficients[2]*x;
    }
    
    /**
     * 
     * @param t1 first point
     * @param incr
     * @param values values at t, t + incr, t + 2*incr
     * @return
     */
    public static double[] coefficients(double t1, double incr, double[] values){
        double[] coefficients = new double[3];
        coefficients[2] = (values[0] + values[2] - 2*values[1])/(2*incr*incr);
        coefficients[1] = (values[2]-values[1])/incr - coefficients[2]*(2*t1+3*incr);
        coefficients[0] = ((t1+2*incr)*values[0]-t1*values[2])/(2*incr) + coefficients[2]*t1*(t1+2*incr);        
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
