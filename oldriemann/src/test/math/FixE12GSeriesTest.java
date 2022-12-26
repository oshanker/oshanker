package math;

import org.junit.Assert;
import org.junit.Test;

public class FixE12GSeriesTest  {
    public void testTestIncrementGValuesAtIndices() {
    }
    
    @Test
    public void testTestChangeToZetaAndDer() {
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
        double[][] ret = fixE12GSeries.testChangeToZetaAndDer(fixE12GSeries.gAtBeta);
        double[] neededZetaIncrement = ret[1];
        double[] actualIncrementInValues = ret[0];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            Assert.assertEquals(actualIncrementInValues[i], neededZetaIncrement[i], 5.0E-6) ;
        }
    }
    
    public void testGradient() {
    }
    
    @Test
    public void testTestChangeToZetaAndDerNoMax() {
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
        double[][] ret = fixE12GSeries.testChangeToZetaAndDerNoMax(fixE12GSeries.gAtBeta);
        double[] neededZetaIncrement = ret[1];
        double[] actualIncrementInValues = ret[0];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            Assert.assertEquals(actualIncrementInValues[i], neededZetaIncrement[i], 5.0E-6) ;
        }
    }
}