package riemann;

public class GseriesSpline extends BaseNormalizedSpline {
    final double[][] y;

    public GseriesSpline(double[][] y) {
        super(y.length);
        this.y = y;
        fit();
    }
    
    @Override
    protected void initSystem(double[] diag, double[] rhs) {
        // TODO Auto-generated method stub

    }

}
