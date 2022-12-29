package math;

import org.junit.Assert;
import org.junit.Test;
import riemann.Interpolate;

import static math.AnalyzeE12GSeries.testGetSavedGSeries1;

public class AnalyzeE12GSeriesTest  {
    @Test
    public void testChangeToZetaAndDer() {
        AnalyzeE12GSeries analyzeE12GSeries = new AnalyzeE12GSeries();
        analyzeE12GSeries.testChangeToZetaAndDer();
    }
    
    @Test
    public void testTestGetSavedGSeries1() {
        double firstZero = 248;
        final double stopValue = 243839.0;
        boolean ignoreMax = true;
        GSeries gAtBeta = Interpolate.readGSeries();
        testGetSavedGSeries1(firstZero, Interpolate.zeroIn, gAtBeta,
            2000002, stopValue, ignoreMax);
    }
}