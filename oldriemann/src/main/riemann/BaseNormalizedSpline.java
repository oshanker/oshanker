package riemann;

public abstract class BaseNormalizedSpline {
    /**
     * si is double the slope.
     */
    final double[] si;
    final int N;
    public BaseNormalizedSpline( int N) {
        super();
        this.si = new double[N-1];
        this.N = N;
    }

    public void fit(){
        double[] diag = new double[N-1];
        double[] rhs = new double[N-1];
        //init
        initSystem(diag, rhs);
        
        si[0] = Double.NEGATIVE_INFINITY;
        
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
    }

    protected abstract void initSystem(double[] diag, double[] rhs);

}
