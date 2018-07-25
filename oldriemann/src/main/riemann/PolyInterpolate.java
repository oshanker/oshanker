package riemann;

public class PolyInterpolate {
    
    public static class BasePolySecondDer extends Interpolate.Poly3{
        double C;
        public BasePolySecondDer(double z0, double z1, double d0, double d1, double secondDer) {
            super(z0, z1, d0, d1);
            C = secondDer/2 + (2*d0+d1)/(z1-z0);
        }

        public double eval2(double x){
            double mult = (x-z0)*(x-z1);
            mult = mult*mult/denom;
            double ret = super.eval1(x)+C*mult;
            return ret;
        }        
        public double der(double x){
            double ret = super.der(x);
            ret += 2*C*(x-z0)*(x-z1)*(2*x-z1-z0)/denom;
            return ret;
        }

        @Override
        void estimateC(double xmin) {
            // dummy
            throw new IllegalStateException("not implemented");
         }
        
    }
    
    public static class Poly5 extends BasePolySecondDer{
        double positionMax;
        double max;

        public Poly5(double z0, double z1, double d0, double d1, 
                double secondDer, double max) {
            super(z0, z1, d0, d1, secondDer);
            if(d0<0){max = -max;}
            this.max= max;
            positionMax = processMax();
        }
        public double eval(double x){
            throw new IllegalStateException("not implemented");
        }        
        public double der(double x){
            throw new IllegalStateException("not implemented");
        }
        
    }

    public static void main(String[] args) {
        BasePolySecondDer basePolySecondDer = new Poly5(1, 2, -1, 1, 0, -0.25); 
        System.out.println(basePolySecondDer.eval1(0));
        System.out.println(basePolySecondDer.eval2(0) );
        System.out.println( " der " +basePolySecondDer.der(1));
        System.out.println(basePolySecondDer.eval2(1.5+Math.sqrt(1.25)));

    }

}
