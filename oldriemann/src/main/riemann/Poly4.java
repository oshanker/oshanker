package riemann;

public class Poly4 extends Poly3 {
   double C;
   double positionMax;

   @Override
   public double getPositionMax() {
      return positionMax;
   }

   double max;
   double cdenom;
   final double sum2;

   public Poly4(double z0, double z1, double d0, double d1, double max) {
      super(z0, z1, d0, d1);
      this.max = max;
      positionMax = processMax();
      cdenom = 4 * C / denom;
      sum2 = (z1 + z0) / 2.0;
   }

   void estimateC(double xmin) {
      double mult = (xmin - z0) * (xmin - z1);
      mult = mult * mult / denom;
      C = (max - eval1(xmin)) / mult;
      cdenom = 4 * C / denom;
   }

   public double eval(double x) {
      double mult = (x - z0) * (x - z1);
      mult = mult * mult / denom;
      double ret = super.eval1(x) + C * mult;
      return ret;
   }

   public double der(double x) {
      double ret = super.der(x);
      //ret += (x-z0)*(x-z1)*(x-sum2)*cdenom;
      ret += 2 * C * (x - z0) * (x - z1) * (2 * x - z1 - z0) / denom;
      return ret;
   }

   public double secondDer(double x) {
      double ret = super.secondDer(x);
      ret += ((x - z0) * (x - z1) + 2 * (x - sum2) * (x - sum2)) * cdenom;
      return ret;
   }

   public double secondDerRHS() {
      double ret = super.secondDerRHS();
      return ret += 2 * C;
   }
}
