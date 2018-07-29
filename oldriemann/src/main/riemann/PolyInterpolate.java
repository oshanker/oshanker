package riemann;

import riemann.Rosser.ZeroInfo;

public class PolyInterpolate {
    
    public static class BasePolySecondDer extends Interpolate.Poly3{
        double C;
        public BasePolySecondDer(double z0, double z1, double d0, double d1, double secondDer) {
            super(z0, z1, d0, d1);
            C = secondDer/2 + (2*d0+d1)/(z1-z0);
        }

        public double eval(double x){
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
        public double secondDerRHS(){
            double ret = super.secondDerRHS();
            return ret += 2*C;        
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
        double D;

        public Poly5(ZeroInfo zeroInput, double secondDer){
            this(zeroInput.lastZero[0],zeroInput.nextValues[0],
                    zeroInput.lastZero[1],zeroInput.nextValues[1],
                    secondDer, zeroInput.lastZero[2]);
            if(!Double.isFinite(D)){
                System.out.println(": " +  z0 +  ", " +  z1 +  ", " +  d0 +  ", " +  d1 +  
                        ", " +  secondDer +  ", " +  max);
                throw new IllegalStateException("NaN"); 
             }
        }

        public Poly5(double z0, double z1, double d0, double d1, 
                double secondDer, double max) {
            super(z0, z1, d0, d1, secondDer);
            if(d0<0){max = -max;}
            this.max= max;
            positionMax = processMax();
            if(positionMax<=z0 || positionMax >= z1){
                System.out.println(": " +  z0 +  ", " +  z1 +  ", " +  d0 +  ", " +  d1 +  
                        ", " +  secondDer +  ", " +  max);
                throw new IllegalStateException("bad data");
            }
        }
        public double eval(double x){
            double mult = (x-z0)*(x-z1);
            mult = mult*mult*(x-z0)/denom;
            double ret = super.eval(x)+D*mult;
            return ret;
        }        
        public double der(double x){
            double ret = super.der(x);
            ret += 2*D*(x-z0)*(x-z0)*(x-z0)*(x-z1)/denom
                    + 3*D*(x-z0)*(x-z0)*(x-z1)*(x-z1)/denom;
            return ret;
        }
        public double secondDerRHS(){
            double ret = super.secondDerRHS();
            return ret += 2*D*h;        
        }

        @Override
        void estimateC(double xmin) {
            double mult = (xmin-z0)*(xmin-z1);
            mult = mult*mult*(xmin-z0)/denom;
            D = (max-super.eval(xmin))/mult;
        }
        
    }

    public static void main(String[] args) {
        double x0 = 2, x1 = 5;
        double N = 7;
        double incr = (x1-x0)/(N-1) ;
//        for (int i = 0; i < N; i++) {
//            double x = x0+incr*i;
//            System.out.println(x + ", " + (Math.pow(x, 5) - 32 - (x-2)*1031));
//            //System.out.println(x + ", der " + (5*Math.pow(x, 4)  - 1031));
//        }
        //(y^5-2^5)-(y-2)*1031
        //5x^4-1031
        double xmin = Math.sqrt(1031.0/5);
        xmin = Math.sqrt(xmin);
        double minValue = Math.pow(xmin, 5) - 32 - (xmin-2)*1031;
        System.out.println(xmin + ", " + minValue);
        System.out.println(xmin + ", der " + (5*Math.pow(xmin, 4)  - 1031));
        Poly5 poly5 = new Poly5(2, 5, -951.0, 2094.0, 160,
                -minValue);
        System.out.println();
        System.out.println(poly5.positionMax +  ", " + poly5.C 
                + ", " + poly5.D + ", secod der " + poly5.secondDerRHS());
        for (int i = 0; i < N; i++) {
            double x = x0+incr*i;
            System.out.println(x + ", " + poly5.eval(x) + ", " + poly5.der(x));
        }
    }

}
