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
        private double D;

        public Poly5(double z0, double z1, double d0, double d1, 
                double secondDer, double max) {
            super(z0, z1, d0, d1, secondDer);
            if(d0<0){max = -max;}
            this.max= max;
            positionMax = processMax();
        }
        public double eval(double x){
            double mult = (x-z0)*(x-z1);
            mult = mult*mult*(x-z0)/denom;
            double ret = super.eval2(x)+D*mult;
            return ret;
        }        
        public double der(double x){
            throw new IllegalStateException("not implemented");
        }
        
    }

    public static void main(String[] args) {
        double x0 = 2, x1 = 5;
        for (int i = 0; i < 7; i++) {
            double x = x0+0.5*i;
            System.out.println(x + ", " + (Math.pow(x, 5) - 32 - (x-2)*1031));
        }
        //4x^4-1031
        double xmin = Math.sqrt(1031.0/4);
        xmin = Math.sqrt(xmin);
        System.out.println(xmin + ", " + (Math.pow(xmin, 5) - 32 - (xmin-2)*1031));
        
        

    }

}
