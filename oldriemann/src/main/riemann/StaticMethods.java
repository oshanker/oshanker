package riemann;

import java.util.function.Function;

public class StaticMethods {
   static void fixLimits(double[] oldest, double[] upper, double xmin, double dermin) {
       if(Math.signum(dermin) == Math.signum(oldest[1])){
           oldest[0]=xmin;
           oldest[1]=dermin;
       } else {
           upper[0]=xmin;
           upper[1]=dermin;
       }
   }

   static void xmin(double[] oldest, double[] upper,
           double[] wts, double precision,
           Function<Double, Double> derivativeFunction){
       double x0 = oldest[0];
       double x1 = upper[0];
       double xmin = (wts[0]*x0 + wts[1]*x1)/(wts[0]+wts[1]);
       double dermin = derivativeFunction.apply(xmin);
       if(Math.abs(dermin)<precision){
           oldest[0]=xmin;
           oldest[1]=dermin;
           oldest[2]=100;
           return;
       }
       fixLimits(oldest, upper, xmin, dermin);
       xmin = (oldest[0] + upper[0])/(2);
       dermin = derivativeFunction.apply(xmin);
       fixLimits(oldest, upper, xmin, dermin);
   }
}
