package math;

import org.junit.Assert;
import org.junit.Test;

import java.util.LinkedList;

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
    
    @Test
    public void testGradient() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
        for (int i = 0; i < FixE12GSeries.desiredSize; i++) {
            zeroInfo.add(fixE12GSeries.nextValues[i]);
        }
        double[] deviation = FixE12GSeries.printZeroInfoWithMax(fixE12GSeries.gAtBeta, zeroInfo);
        Assert.assertTrue(deviation[0] < 1.0E-9);
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