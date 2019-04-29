package riemann;

/**
 * Not used
 * @author shankero
 *
 */
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
