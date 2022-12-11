package riemann;

public class Poly7 {
    final double a, b, c;
    double d0, d1, d2;
    
    public Poly7(double a, double b, double c,
                 double a1, double b1, double c1) {
        this.a = a;
        this.b = b;
        this.c = c;
        d0 = a1;
        d1 = b1;
        d2 = c1;
    }
    
    double der(double x) {
        if (Math.abs(x-a) < 1.0E-10 ) {
            return d0;
        }
        if (Math.abs(x-b) < 1.0E-10 ) {
            return d1;
        }
        if ( Math.abs(x-c) < 1.0E-10) {
            return d2;
        }
        double prod = (x-a)*(x-b)*(x-c);
        //double ret = (x-a)*(x-b) + (x-b)*(x-c) + (x-a)*(x-c);
        double ret = prod* (1/(x-a) + 1/(x-b) + 1/(x-c));
        return ret;
    }
    
    double eval(double x) {
        double prod = (x-a)*(x-b)*(x-c);
        double ret = prod;
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
                return x0;
            }
        }
        return next;
    }
    
    public static void main(String[] args) {
        Poly7 poly7term = new Poly7(0, 1, 2, 2, -1, 2);
        tabulate(poly7term);
        poly7term.positionMax(0.5, 0, 1);
    }
    
    private static void tabulate(Poly7 poly7term) {
        for (double x = -0.1; x < 2.2; x += 0.05) {
            System.out.println(x +
                " " + poly7term.eval(x) +
                " " + poly7term.der(x)
            );
        }
    }
}
