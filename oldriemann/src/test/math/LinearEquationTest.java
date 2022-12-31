package math;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.LinkedList;

import static riemann.StaticMethods.changeToDer;
import static riemann.StaticMethods.changeToZeta;
import static riemann.StaticMethods.changeToZetaAndDer;
import static riemann.StaticMethods.evaluateAtT;
import static org.junit.Assert.assertEquals;
import static riemann.StaticMethods.findFile;
import static riemann.StaticMethods.getSavedGSeries;
import static riemann.StaticMethods.gramE12;

public class LinearEquationTest {
    static final int initialPadding = 40;
    static double[] nextValues = {
        243831.456494008, -22.69554476177354, 1.538114456203189};
    static double[] nextValues1 = {243831.660443468, 28.68660904273845, 3.261849766062147};
    static double[] nextValues2 = {243831.92506103282, -46.745064213360436, 4.265426650034286};
    
    static double[] pointBeingInflunced = {nextValues[0], nextValues1[0], nextValues2[0]};
    static GSeries gAtBeta = null;
    static double[] initial = null;
    int midIdxCausingInfluence = 9994;
    
    @BeforeClass
    public static void beforeClass() throws Exception {
        gAtBeta = getgSeries(pointBeingInflunced[0]);
        initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        //System.out.println("initial " + Arrays.toString(initial));
    }
    
    @Test
    public void testIncrementGValuesAtIndices() {
    
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
    
        // row = gseries indices NO!
        // col = points being influenced
        double[][] coefficients = new double[][]{zetaCoeff[0], derCoeff[0]};
        LinearEquation linearEquation = new LinearEquation(coefficients);
    
        double[] solution = linearEquation.solve(new double[]{1, 0.5});
        System.out.println("Required g increment " + Arrays.toString(solution));
    
        System.out.println("multiply(coefficients, solution) " +
            Arrays.toString(LinearEquation.multiply(coefficients, solution)));
        double[][] coefficients1 = new double[][]{zetaCoeff[1], derCoeff[1]};
        System.out.println("multiply(coefficients1, solution) " +
            Arrays.toString(LinearEquation.multiply(coefficients1, solution)));
        double[][] requiredGIncrements = {
            {0.13861482493261557, 0.8115131385775348},
            {0, 0}
        };
        gAtBeta.incrementGValuesAtIndices(midIdxCausingInfluence, requiredGIncrements);
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("after " + Arrays.toString(after));
        double[] actualIncrement = new double[after.length];
        for (int i = 0; i < actualIncrement.length; i++) {
            actualIncrement[i] = after[i] - initial[i];
        }
        assertEquals(1.0, actualIncrement[0], 0.000001);
        assertEquals(0.5, actualIncrement[1], 0.000001);
    
        gAtBeta.decrementGValuesAtIndices(midIdxCausingInfluence, requiredGIncrements);
        System.out.println("zetaCoeff " + Arrays.deepToString(zetaCoeff));
        System.out.println("derCoeff " + Arrays.deepToString(derCoeff));
    
        //================================
        int[] indices = {midIdxCausingInfluence, midIdxCausingInfluence + 1};
        double[][] zetaDerCoeff = changeToZetaAndDer(
            gAtBeta,
            initialPadding,
            pointBeingInflunced,
            initial,
            indices,
            0.125
        );
        System.out.println("zetaderCoeff ");
        LinearEquation.printMatrix(zetaDerCoeff);
        double[] testGIncrements = {
            0.13861482493261557, 0.8115131385775348,
            0, 0
        };
        System.out.println("test prod ");
        System.out.println("actualIncrement " + Arrays.toString(actualIncrement));
        double[] expected = LinearEquation.multiply(zetaDerCoeff, testGIncrements);
        System.out.println(Arrays.toString(expected));
        for (int i = 0; i < solution.length; i++) {
            Assert.assertEquals(expected[i], actualIncrement[i], 1.0E-9);
        }
    }
    
    @Test
    public void testoverdetermined() {
        double[][] a = {
            {1, 1},
            {1, 2},
            {1, 3},
        };
        double[][] at = LinearEquation.transpose(a);
        double[] y = {6, 11, 16};
        double[][] coefficients = LinearEquation.multiply(at, a);
        double[] values = LinearEquation.multiply(at, y);
        LinearEquation linearEquation = new LinearEquation(coefficients);
        double[] solution = linearEquation.solve(values);
        System.out.println(Arrays.toString(solution));
        double[] expected = {1, 5};
        for (int i = 0; i < solution.length; i++) {
            Assert.assertEquals(expected[i], solution[i], 1.0E-9);
            
        }
    }
    
    @Test
    public void testCoeff1() {
        double[][] coefficients = {
            {1, -1},
            {1, 1},
        };
        LinearEquation linearEquation = new LinearEquation(coefficients);
        LinearEquation.printMatrix(linearEquation.coefficients);
        Assert.assertEquals(2.0, linearEquation.determinant(), 1.0E-9);
        System.out.println("=======");
        coefficients = new double[][]{
            {1, 0.999},
            {1, 1},
        };
        linearEquation = new LinearEquation(coefficients);
        LinearEquation.printMatrix(linearEquation.coefficients);
        Assert.assertEquals(0.001, linearEquation.determinant(), 1.0E-9);
    }
    
    @Test
    public void testCoeff() {
        int n = 5;
        double[][] coefficients = new double[n][n];
        double[] values = new double[n];
        for(int row = 0; row < n; row++)
        {
            coefficients[row][row ] = row + 1;
            for(int col = row+1; col < n; col ++)
            {
                coefficients[row][col ] = 1;
            }
            values[row] = 2*row + 1;
        }
        double[][] transpose = LinearEquation.transpose(coefficients);
        LinearEquation linearEquation = new LinearEquation(transpose);
        LinearEquation.printMatrix(transpose);
        System.out.println("=======");
        LinearEquation.printMatrix(linearEquation.coefficients);
        double[] solution = linearEquation.solve(values);
        System.out.println(Arrays.toString(solution));
        
    }
    
    @Test
    public void testChangeToZetaAndDer() {
        int[] indices = {midIdxCausingInfluence, midIdxCausingInfluence+1, midIdxCausingInfluence+2};
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
        LinearEquation linearEquation = new LinearEquation(zetaDerCoeff );
    
        double[] solution = linearEquation.solve(
            new double[]{1, 0.5, -0.16821753315955723, -2.0541382955436056, -0.06545828959341304, -1.0119912921417438});
        System.out.println("Required g increment " );
        System.out.println( Arrays.toString(solution));
        double[] testGIncrements = {
            0.13861482493261557, 0.8115131385775348,
            0, 0, 0, 0
        };
        for (int i = 0; i < solution.length; i++) {
            Assert.assertEquals(testGIncrements[i], solution[i], 1.0E-9);
        
        }
    }
    
    private static GSeries getgSeries(double pointBeingInflunced) throws IOException {
        int fileIdx = findFile(pointBeingInflunced);

        double t0 = gramE12[fileIdx][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = getSavedGSeries(t0, offset);
        return gAtBeta;
    }


}
