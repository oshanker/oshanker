package riemann;

public class Spline {
    

    public static class Splinei {
        public final double xi, h;
        public final double yi, yi1;
        public final double si, si1;
        //first argument should be 0, (zero based)!!!
        // or maybe, just take h as the increment?
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

        public double seconddDer(double x){
            double ret = 0;
            double xifac = (x-xi)/h;
            ret = 4*xifac*((si1-2*yi1/h))/h;
            ret += 2*(yi1 + (x-xi-h)*(si1-2*yi1/h))/(h*h);
            double xi1fac = (x-xi-h)/h;
            ret += 4*xi1fac*((si+2*yi/h))/h;
            ret += 2*(yi + (x-xi)*(si+2*yi/h))/(h*h);
            return ret;
        }

    }

    public static void main(String[] args) {
        double xi = 0, h = 8;
        double yi = 0, yi1 = 0;
        double slopei = 0, slopei1 = 64;
        //x^3-8x^2
        Splinei splinei = new Splinei(xi, h, yi, yi1, slopei, slopei1);
        int N = 5;
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
        System.out.println("From second der at x = 4 :" + (2*(si[3]+4*si[2]+si[1])) + ", " + (3*(y[3]-y[1])));
    }

}
