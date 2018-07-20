package riemann;

public class Spline {
    

    public static class Splinei {
        public final double xi, h;
        public final double yi, yi1;
        public final double si, si1;
        public Splinei(double xi, double h, double yi, double yi1,
                 double si, double si1) {
            this.xi = xi;
            this.h = h;
            this.yi = yi;
            this.yi1 = yi1;
            this.si = si;
            this.si1 = si1;
        }
        
        public double eval(double x){
            double ret = 0;
            double xifac = (x-xi)/h;
            ret = xifac*xifac*(yi1 + (x-xi-h)*(si1-2*yi1/h));
            double xi1fac = (x-xi-h)/h;
            ret += xi1fac*xi1fac*(yi + (x-xi)*(si+2*yi/h));
            return ret;
        }
        
        public double der(double x){
            double ret = 0;
            double xifac = (x-xi)/h;
            ret = 2*xifac*(yi1 + (x-xi-h)*(si1-2*yi1/h))/h;
            ret += xifac*xifac*((si1-2*yi1/h));
            double xi1fac = (x-xi-h)/h;
            ret += 2*xi1fac*(yi + (x-xi)*(si+2*yi/h))/h;
            ret += xi1fac*xi1fac*((si+2*yi/h));
            return ret;
        }

    }

    public static void main(String[] args) {
        double xi = 1, h = 1;
        double yi = 2, yi1 = 2;
        double si = -2, si1 = 2;
        //2*(x-1)*(x-2) + 2
        Splinei splinei = new Splinei(xi, h, yi, yi1, si, si1);
        for (int i = 0; i < 3; i++) {
            double x = 1+0.5*i;
            System.out.println(x + ", " + splinei.eval(x) + ", " + splinei.der(x));
        }
    }

}
