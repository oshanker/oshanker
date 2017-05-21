/**
 * 
 */
package math;

/**
 * @author shankero
 *
 */
public class ZeroPoly {
    double[] roots;
    double[] slopes;
    

    public ZeroPoly(double[] roots, double[] slopes) {
        super();
        this.roots = roots;
        this.slopes = slopes;
    }
    
    public double prodExclude(double t, int excludeIndex){
        double val = 1;
        for (int i = 0; i < roots.length; i++) {
            if(i == excludeIndex){continue;}
            val*=(t-roots[i]);
        }
        return val;
    }
    
    public double prod(double t){
        double val = 1;
        for (int i = 0; i < roots.length; i++) {
            val*=(t-roots[i]);
        }
        return val;
    }
    
    public  double[] coefficients(double t1, double incr){
        double[] z = new double[3];
        for (int i = 0; i < z.length; i++) {
            double x  = t1+i*incr;
            z[i] = eval(x);
        }
        double[] coeff = Quadratic.coefficients(t1, incr, z);
        return coeff;
    }   
    
    public double eval(double t){
        double val = 0;
        double prod = prod(t);
        for (int i = 0; i < roots.length; i++) {
            double norm = prodExclude(roots[i], i);
            //System.out.println(i + " norm "+ norm);
            double valTerm = prod*(slopes[i])/((t-roots[i])*norm*norm);
            val += valTerm;
            //System.out.println(valTerm);
        }
        val *= prod;
        return val;
    }


    /**
     * @param args
     */
    public static void main(String[] args) {
        double[] roots = new double[]{0.4375, 0.5061};
        double[] slopes = new double[]{-26.17, 14.50};
        ZeroPoly zeroPoly = new ZeroPoly(roots, slopes);
        double eval = zeroPoly.eval(0.591);
        System.out.println("eval " + eval);
    }

}
