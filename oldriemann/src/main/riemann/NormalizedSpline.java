package riemann;

import java.util.Arrays;

public class NormalizedSpline {
    public static class Splinei {
        public final double xi, originalH;
        public final double yi, yi1;
        public final double si, si1;
        public Splinei(double xi, double h, double yi, double yi1,
                 double si, double si1) {
            this.xi = xi/h;
            this.yi = yi;
            this.yi1 = yi1;
            this.si = si;
            this.si1 = si1;
            this.originalH = h;
        }
        
        public double eval(double x){
            double ret = 0;
            double xifac = (x/originalH-xi);
            double xi1fac = (xifac-1);
            ret = xifac*xifac*(yi1 + (xi1fac)*(si1-2*yi1));
            ret += xi1fac*xi1fac*(yi + xifac*(si+2*yi));
            return ret;
        }
        
        public double der(double x){
            double ret = 0;
            double xifac = (x/originalH-xi);
            double xi1fac = (xifac-1);
            ret = 2*xifac*(yi1 + (xi1fac)*(si1-2*yi1));
            ret += xifac*xifac*((si1-2*yi1));
            ret += 2*xi1fac*(yi + (xifac)*(si+2*yi));
            ret += xi1fac*xi1fac*((si+2*yi));
            return ret;
        }

        public double seconddDer(double x){
            double ret = 0;
            double xifac = (x/originalH-xi);
            double xi1fac = (xifac-1);
            ret = 4*xifac*((si1-2*yi1));
            ret += 2*(yi1 + (xi1fac)*(si1-2*yi1));
            ret += 4*xi1fac*((si+2*yi));
            ret += 2*(yi + (xifac)*(si+2*yi));
            return ret;
        }
    }
    
    final double[] y;
    final int N;

    public NormalizedSpline(double[] y) {
        super();
        this.y = y;
        this.N = this.y.length;
        /**
0, 0.0, 0.0, -1024.0
1, -24.0, -160.0, -256.0
2, -64.0, -128.0, 512.0
3, -72.0, 96.0, 1280.0
4, 0.0, 512.0, 2048.0
5, 200.0, 1120.0, 2816.0
From second der at i = 2 :-115.2, -144.0
From second der at i = 3 :153.6, 192.0
From second der at i = 4 :652.8, 816.0
From first spline :-224.0, -224.0
From last spline :560.0, 560.0
         */
    }
    
    public void fit(){
        //si[1]
        System.out.println();
        double[] si = new double[N-1];
        double[] diag = new double[N-1];
        double[] rhs = new double[N-1];
        //init
        diag[1] = 4.0;
        rhs[1] = (5*y[2]-4*y[1]-y[0]);
        for (int i = 2; i < diag.length-1; i++) {
            diag[i] = 4.0;
            rhs[i] = (3*(y[i+1]-y[i-1]));
        }
        diag[diag.length-1] = 4.0;
        rhs[diag.length-1] = (-5*y[N-3]+4*y[N-2]+y[N-1]);
        
        int index = N-3;
        double mult = 1/diag[index+1];
        diag[index] -= 2.0*mult;
        rhs[index] -= rhs[index+1]*mult;
        index--;
        while(index>=2){
            mult = 1/diag[index+1];
            diag[index] -= mult;
            rhs[index] -= rhs[index+1]*mult;
            index--;
        }
        mult = 2.0/diag[index+1];
        diag[index] -= mult;
        rhs[index] -= rhs[index+1]*mult;
        index--;
        
        index = 1;
        si[index] = rhs[index]/diag[index];
        index++;
        while(index<N-2){
            si[index] = (rhs[index]-si[index-1])/diag[index];
            index++;
        }
        si[index] = (rhs[index]-2*si[index-1])/diag[index];
        
        System.out.println(Arrays.toString(diag));
        System.out.println(Arrays.toString(rhs));
        System.out.println(Arrays.toString(si));
    }


    public static void main(String[] args) {
        //x^3-8x^2
        int N = 6;
        double[] y = new double[N];
        for (int i = 0; i < N; i++) {
            double x = 2*i;
             y[i] = x*x*x - 8*x*x;
            System.out.println(i + ", " + y[i] );
        }
        NormalizedSpline normalizedSpline = new NormalizedSpline(y);
        normalizedSpline.fit();
    }


}
