package math;

import org.junit.Assert;
import org.junit.Test;
import riemann.Interpolate;

import java.util.Arrays;
import java.util.LinkedList;

import static math.FixE12GSeries.TEST_VALES;
import static riemann.StaticMethods.changeToDer;
import static riemann.StaticMethods.changeToZeta;
import static riemann.StaticMethods.changeToZetaAndDer;
import static riemann.StaticMethods.evaluateAtT;

public class FixE12GSeriesTest  {
    static final int initialPadding = 40;
    static double[][] badSigmum = {
        {490.8522279298088, 27.694175364012526, 2.2591989985805814, 490.9649293406849},
        {491.077630751561, -38.40311109971628, -4.196432798508113, 491.20420495587484},
        {491.33077916018874, 72.25062982634255, 14.450288488557158, 0.0},
        {491.89945544848615, -6.8841110643202965, -0.09212140202372199},
    };
    
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
    public void testGradient1() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        for (int i = 0; i < badSigmum.length; i++) {
            zeroInfo.add(badSigmum[i]);
        }
        GSeries gAtBeta = Interpolate.readGSeries();
    
        FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, gAtBeta);
        double[][] ret = FixE12GSeries.printZeroInfoWithMax(fixE12GSeries.gAtBeta, zeroInfo);
        double[] deviation = ret[0];
        Assert.assertTrue(deviation[0] < 1.0E-9);
        Assert.assertEquals(1999906, (int)ret[1][0] );
    }
    
    @Test
    public void testChangeToZetaAndDerNoMax2() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        for (int i = 0; i < badSigmum.length; i++) {
            zeroInfo.add(badSigmum[i]);
        }
        GSeries gAtBeta = Interpolate.readGSeries();
        
        FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, gAtBeta);
        double[][] ret = fixE12GSeries.testChangeToZetaAndDerNoMax(
            fixE12GSeries.gAtBeta, 2024, false);
        double[] deviation = ret[2];
        System.out.println(Arrays.toString(deviation));
        Assert.assertTrue(deviation[0] < 1.0E-9);
    }
    
    @Test
    public void testGradient() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
        for (int i = 0; i < FixE12GSeries.desiredSize; i++) {
            zeroInfo.add(fixE12GSeries.nextValues[i]);
        }
        double[][] ret = FixE12GSeries.printZeroInfoWithMax(fixE12GSeries.gAtBeta, zeroInfo);
        double[] deviation = ret[0];
        Assert.assertTrue(deviation[0] < 1.0E-9);
        Assert.assertEquals(1999906, (int)ret[1][0] );
    }
    
    @Test
    public void testTestChangeToZetaAndDerNoMax1() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        for (int i = 0; i < TEST_VALES.length; i++) {
            zeroInfo.add(TEST_VALES[i]);
        }
    
        GSeries gAtBeta = Interpolate.readGSeries();
    
        FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, gAtBeta);
        double[][] ret = fixE12GSeries.testChangeToZetaAndDerNoMax(
            gAtBeta, -1, true);
        double[] neededZetaIncrement = ret[1];
        double[] actualIncrementInValues = ret[0];
        double[] deviation = ret[2];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            Assert.assertEquals(actualIncrementInValues[i], neededZetaIncrement[i], 5.0E-6) ;
        }
        Assert.assertTrue(deviation[0] < 1.0E-9);
        //System.out.println(" " + Arrays.toString(actualIncrementInValues));
    }
    
    @Test
    public void testTestChangeToZetaAndDerNoMax() {
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
        double[][] ret = fixE12GSeries.testChangeToZetaAndDerNoMax(
            fixE12GSeries.gAtBeta, 1999906, true);
        double[] neededZetaIncrement = ret[1];
        double[] actualIncrementInValues = ret[0];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            Assert.assertEquals(actualIncrementInValues[i], neededZetaIncrement[i], 5.0E-6) ;
        }
        //System.out.println(" " + Arrays.toString(actualIncrementInValues));
    }
    
    @Test
    public void testIncrementGValuesAtIndices() {
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
    
        double[][] nextValues = fixE12GSeries.nextValues;
        double[] pointBeingInflunced = new double[]{
            nextValues[0][0], nextValues[1][0],
            nextValues[2][0], nextValues[3][0],
            nextValues[4][0],
        };
        GSeries gAtBeta = fixE12GSeries.gAtBeta;
        double[] initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("initial " + Arrays.toString(initial));
    
        int midIdxCausingInfluence = fixE12GSeries.midIdxCausingInfluence;
        double[][] zetaCoeff = changeToZeta(
            gAtBeta,
            initialPadding,
            pointBeingInflunced,
            initial,
            midIdxCausingInfluence,
            0.125
        );
        double[][] derCoeff = changeToDer(
            gAtBeta, initialPadding, pointBeingInflunced, initial,
            midIdxCausingInfluence, 0.125
        );
        
        // row = gseries indices
        // col = points being influenced
        double[][] coefficients = new double[][]{zetaCoeff[2], derCoeff[2]};
        LinearEquation linearEquation = new LinearEquation(coefficients );
        
        double[] solution = linearEquation.solve(new double[]{1, 0.5});
        System.out.println("Required g increment " + Arrays.toString(solution));
        
        System.out.println("multiply(coefficients, solution) " +
            Arrays.toString(LinearEquation.multiply(coefficients, solution)));
        double[][] coefficients1 = new double[][]{zetaCoeff[1], derCoeff[1]};
        System.out.println("multiply(coefficients1, solution) " +
            Arrays.toString(LinearEquation.multiply(coefficients1, solution)));
        double[][] requiredGIncrements = {
            {
                -2.800779362553829, -0.7056081953041438,
                //-0.6298688316508069, 0.2743224605076784,
            },
            {0, 0}
        };
        gAtBeta.incrementGValuesAtIndices(midIdxCausingInfluence, requiredGIncrements);
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("after " + Arrays.toString(after));
        double[] actualIncrement = new double[after.length];
        for (int i = 0; i < actualIncrement.length; i++) {
            actualIncrement[i] = after[i] - initial[i];
        }
        System.out.println("actualIncrement " );
        System.out.println( Arrays.toString(actualIncrement));
        
        gAtBeta.decrementGValuesAtIndices(midIdxCausingInfluence, requiredGIncrements);
        System.out.println("zetaCoeff " + Arrays.deepToString(zetaCoeff));
        System.out.println("derCoeff " + Arrays.deepToString(derCoeff));
        
        //================================
        int[] indices = {midIdxCausingInfluence, midIdxCausingInfluence+1};
        double[][] zetaDerCoeff = changeToZetaAndDer(
            gAtBeta,
            initialPadding,
            pointBeingInflunced,
            initial,
            indices,
            0.125
        );
        System.out.println("zetaderCoeff " );
        LinearEquation.printMatrix(zetaDerCoeff);
        Assert.assertEquals(1.0, actualIncrement[4], 1.0E-6);
        Assert.assertEquals(0.5, actualIncrement[5], 1.0E-6);
    }
    
}