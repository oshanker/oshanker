package math;

import org.junit.Assert;
import org.junit.Test;
import riemann.Interpolate;

import java.util.Arrays;
import java.util.LinkedList;

import static math.FixE12GSeries.TEST_VALES;

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
    public void testTestChangeToZetaAndDerNoMax1() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        for (int i = 0; i < TEST_VALES.length; i++) {
            zeroInfo.add(TEST_VALES[i]);
        }
    
        GSeries gAtBeta = Interpolate.readGSeries();
    
        FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, -1000, gAtBeta);
        double[][] ret = fixE12GSeries.testChangeToZetaAndDerNoMax(gAtBeta);
        double[] neededZetaIncrement = ret[1];
        double[] actualIncrementInValues = ret[0];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            Assert.assertEquals(actualIncrementInValues[i], neededZetaIncrement[i], 5.0E-6) ;
        }
        //System.out.println(" " + Arrays.toString(actualIncrementInValues));
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
        //System.out.println(" " + Arrays.toString(actualIncrementInValues));
    }
}