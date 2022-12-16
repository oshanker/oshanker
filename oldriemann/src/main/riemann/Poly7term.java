package riemann;

public class Poly7term {
    final double a, b, c;
    double A, B;
    
    public Poly7term(double a, double b, double c, double a1, double b1) {
        this.a = a;
        this.b = b;
        this.c = c;
        A = a1;
        B = b1;
    }
    
    double der(double x) {
        if (Math.abs(x-a) < 1.0E-10 || Math.abs(x-b) < 1.0E-10 | Math.abs(x-c) < 1.0E-10) {
            return 0;
        }
        double prod = (x-a)*(x-b)*(x-c);
        prod = prod*prod;
        double ret = 2*(A*x+B)*prod*
            (1/(x-a) + 1/(x-b) + 1/(x-c)) + A*prod;
        return ret;
    }
    
    double eval(double x) {
        double prod = (x-a)*(x-b)*(x-c);
        prod = prod*prod;
        double ret = prod*(A*x+B);
        return ret;
    }
    
    double positionMax(double x0, double xa, double xb) {
        double der = der(x0);
        System.out.println(x0 + " " + der);
        if (Math.abs(der) < 1.0E-8) {
            return x0;
        }
        double incr = 0.001*(xb-xa);
        double next = x0 + incr;
        double derNext = der(next);
    
        for (int i = 0; i < 4; i++) {
            double oldx0 = x0;
            double oldder = der;
            der = derNext;
            x0 = next;
            next = oldx0 - oldder*(next-oldx0)/(derNext-oldder);
            derNext = der(next);
            System.out.println(next + " " + derNext);
            if (Math.abs(derNext) < 1.0E-8) {
                return next;
            }
        }
        return next;
    }
    
    public static void main(String[] args) {
        Poly7term poly7term = new Poly7term(0, 1, 2, 1, 0);
        tabulate(poly7term);
        poly7term.positionMax(0.5, 1, 0);
    }
    
    private static void tabulate(Poly7term poly7term) {
        for (double x = -0.1; x < 2.2; x += 0.05) {
            System.out.println(x +
                " " + poly7term.eval(x) +
                " " + poly7term.der(x)
            );
        }
    }
}
