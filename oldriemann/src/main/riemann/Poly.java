package riemann;

public interface Poly {
    double der(double x);
    
    double eval(double x);
    
    double getPositionMax();
    
    double secondDer(double upperLimit);
}
