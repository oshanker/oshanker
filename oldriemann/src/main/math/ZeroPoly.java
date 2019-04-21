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
    
    public double prodExclude(double t, int excludeIndex0, int excludeIndex1){
        if(roots.length == 2) {return 0;}
    	double val = 1;
        for (int i = 0; i < roots.length; i++) {
            if(i == excludeIndex0 || i == excludeIndex1){continue;}
            val*=(t-roots[i]);
        }
        return val;
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
    
    public double secondDer(int index) {
    	double ret = 0;
    	double t = roots[index];
        for (int j = 0; j < roots.length; j++) {
            if(j == index){continue;}
            ret += 1/(t-roots[j]);
        }
    	//System.out.println(ret);
    	ret *= 2*slopes[index];
    	//System.out.println(ret);
    	double peindex = prodExclude(t, index);
        for (int j = 0; j < roots.length; j++) {
            if(j == index){continue;}
        	double pej = prodExclude(roots[j], j);
            //System.out.println(j + " norm "+ pej);
        	ret += slopes[j]*peindex*(prodExclude(t, index, j))/(pej*pej);
        }
        return 2*ret;
    }
    
    public double der(int index){
        double val = 0;
    	double t = roots[index];
        double prod = prodExclude(t, index);
        for (int i = 0; i < roots.length; i++) {
        	//only one term contributes
            double norm = prodExclude(roots[i], i);
            //System.out.println(i + " norm "+ norm);
            double valTerm = prodExclude(t, i)*(slopes[i])/(norm*norm);
            val += valTerm;
            //System.out.println(i + ", " + valTerm);
        }
        val *= prod;
        return val;
    }
    
    public double eval(double t){
        double val = 0;
        double prod = prod(t);
        for (int i = 0; i < roots.length; i++) {
            double norm = prodExclude(roots[i], i);
            //System.out.println(i + " norm "+ norm);
            double valTerm = (slopes[i])/((t-roots[i])*norm*norm);
            val += valTerm;
            //System.out.println(valTerm);
        }
        val *= prod*prod;
        return val;
    }


    /**
     * @param args
     */
    public static void main(String[] args) {
        double[] roots = new double[]{1, 2, 3};
        double[] slopes = new double[]{2, -1, 2};
        ZeroPoly zeroPoly = new ZeroPoly(roots, slopes);
        double eval = zeroPoly.eval(0.0);
        System.out.println("eval " + eval);
        for (int i = 0; i < slopes.length; i++) {
            double der = zeroPoly.secondDer(i);
            System.out.println( " der " + der);
		}
    }

}
