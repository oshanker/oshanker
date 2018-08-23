package riemann;

public class TwoPointFit {
    private final double x0,  x1;
    private final double A;
    private final double B;
    final double C;
    final double D;
    final double E;
    final double F;
    public TwoPointFit(double x0, double y0, double d0, double dd0,
            double x1, double y1, double d1, double dd1)
    {
        this.x0 = x0;
        this.x1 = x1;
        double h = x1-x0;
        double h2 = h*h;
        A = y1/h;
        B = -y0/h;
        C = (d1-A-B)/h2;
        D = (d0-A-B)/h2;
        double h3 = h2*h;
        E = (dd1-2*h*(2*C+D))/(2*h3);
        F = (dd0+2*h*(C+2*D))/(-2*h3);
        System.out.println("A " + A + " B " + B);
        System.out.println("C " + C + " D " + D);
        System.out.println("E " + E + " F " + F);
    }
    
    public double eval(double x){
        double f0 = (x-x0), f1 = (x-x1);
        double ret = A*f0 + B*f1;
        double prod = f0*f1;
        ret += prod*(C*f0 + D*f1);
        ret += prod*prod*(E*f0 + F*f1);
        return ret;
    }
    
    /*
h = x1-x0
 y1*(x-x0)/h - y0*(x-x1)/h //val
+ (x-x0)(x-x1)*( (x-x0)(d1+y0/h-y1/h)/h^2 + (x-x1)(d0-y1/h+y0/h)/h^2 )//der
+ (x-x0)(x-x1)*(x-x0)(x-x1)*( E(x-x0) + F(x-x1) ) //second der
LHS -h*2*C +xxx 2*(-h)*D; 2*(-h^3)*F
RHS 2*h*C + 2*h*D; 2*(h^3)*E

     */

    public static void main(String[] args) {
        // x^3-8 - (x-2)*(125-8)/3 + x^4
        // 3x^2-39 + 4*x^3
        // 6x + 12*x^2
        TwoPointFit twoPointFit = new TwoPointFit(0, 70, -39, 0, 
                2, 0 + 16, -27 + 32, 12 + 48);
        System.out.println(twoPointFit.eval(5));//625
        System.out.println(twoPointFit.eval(1));//33
    }

}
