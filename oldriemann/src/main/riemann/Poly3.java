package riemann;

public abstract class Poly3 {
   final double z0, z1;
   final double d0, d1;
   double denom;
   final double h;

   final double pow2;
   final double pow1;
   final double pow0;
   final double cross;

   Poly3(double z0, double z1, double d0, double d1) {
      this.z0 = z0;
      this.z1 = z1;
      this.d0 = d0;
      this.d1 = d1;
      h = (z1 - z0);
      denom = h * h;

      pow2 = (d1 + d0) / denom;
      pow0 = ((d1 + d0) * (2 * z0 * z1) + z0 * z0 * d1 + z1 * z1 * d0) / denom;
      cross = (z0 * d1 + z1 * d0) / denom;
      pow1 = 2 * (pow2 * (-z0 - z1) - cross);
   }

   public double eval1(double x) {
      double ret = (x - z0) * (x - z1) / denom;
      ret *= ((x - z0) * d1 + (x - z1) * d0);
//            double ret = (x-z0)*(x-z1);
//            ret *= pow2*x - cross;
      return ret;
   }

   public double der(double x) {
      double ret = (
            (d1 + d0) * (x - z0) * (x - z1) + (x - z0) * ((x - z0) * d1
                  + (x - z1) * d0) + (x - z1) * ((x - z0) * d1 + (x - z1) * d0)
      ) / denom;
//            double ret = (3*x*x*pow2
//                    + x*pow1
//                    + pow0);
      return ret;
   }

   public double secondDer(double x) {
      double ret = (6 * x * pow2
            + pow1);
      return ret;
   }

   public double secondDerRHS() {
      double ret = 2 * (d0 + 2 * d1) / h;
      return ret;
   }

   abstract void estimateC(double xmin);

   public abstract double eval(double x);

   public abstract double getPositionMax();

   protected double processMax() {
      double[] oldest = new double[]{z0, d0, 0};
      double[] upper = new double[]{z1, d1, 1};
      double[] wts = new double[]{Math.abs(d1), Math.abs(d0)};
      double xmin = z0 - d0 * (z1 - z0) / (d1 - d0);
      double precision = 0.001 * Math.abs(d0);
      for (int i = 0; i < 10; i++) {
         StaticMethods.xmin(oldest, upper, wts, precision, (x) -> der(x));
         if (oldest[2] > 99) {
            xmin = oldest[0];
            break;
         }
         wts[0] = Math.abs(upper[1]);
         wts[1] = Math.abs(oldest[1]);
         if (wts[0] > wts[1]) {
            xmin = oldest[0];
         } else {
            xmin = upper[0];
         }
      }
      estimateC(xmin);
      oldest = new double[]{z0, d0, 0};
      upper = new double[]{z1, d1, 1};
      double dermin = der(xmin);
      StaticMethods.fixLimits(oldest, upper, xmin, dermin);
      wts = new double[]{Math.abs(upper[1]), Math.abs(oldest[1])};
      precision = 0.0000001;
      for (int i = 0; i < 10; i++) {
         estimateC(xmin);
         StaticMethods.xmin(oldest, upper, wts, precision, (x) -> der(x));
         if (oldest[2] > 99) {
            xmin = oldest[0];
            break;
         }
         wts[0] = Math.abs(upper[1]);
         wts[1] = Math.abs(oldest[1]);
         if (wts[0] > wts[1]) {
            xmin = oldest[0];
         } else {
            xmin = upper[0];
         }
      }
      return xmin;
   }

}
