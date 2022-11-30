package math;

import org.junit.Test;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;

import static math.FixGSeries.changeToDer;
import static math.FixGSeries.changeToZeta;
import static riemann.StaticMethods.evaluateAtT;
import static org.junit.Assert.assertEquals;
import static riemann.StaticMethods.findFile;
import static riemann.StaticMethods.getSavedGSeries;
import static riemann.StaticMethods.gramE12;

public class LinearEquationTest {
    final int initialPadding = 40;

    @Test
    public void testIncrementGValuesAtIndices() throws IOException {
        double[] nextValues = {
                243831.456494008, -22.69554476177354, 1.538114456203189};
        double[] nextValues1 = {243831.660443468, 28.68660904273845, 3.261849766062147} ;
        double[] pointBeingInflunced = {nextValues[0], nextValues1[0]};
        GSeries gAtBeta = getgSeries(pointBeingInflunced[0]);
        double[] initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("initial " + Arrays.toString(initial));
        int midIdxCausingInfluence  = 9994;

        double[] zetaCoeff = changeToZeta(
                gAtBeta,
                initialPadding,
                pointBeingInflunced[0],
                initial[0],
                midIdxCausingInfluence,
                0.125
        );
        System.out.println("zetaCoeff " + Arrays.toString(zetaCoeff));
        double[] derCoeff = changeToDer(
                gAtBeta, initialPadding, pointBeingInflunced[0], initial[1],
                midIdxCausingInfluence, 0.125
        );
        System.out.println("derCoeff " + Arrays.toString(derCoeff));

        // row = gseries indices
        // col = points being influenced
        double[][] coefficients = new double[][]{zetaCoeff,derCoeff};
        LinearEquation linearEquation = new LinearEquation(coefficients );

        double[] solution = linearEquation.solve(new double[]{1, 0.5});
        System.out.println("Required g increment " + Arrays.toString(solution));

        System.out.println("multiply(coefficients, solution) " + Arrays.toString(LinearEquation.multiply(coefficients, solution)));
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
        System.out.println("actualIncrement " + Arrays.toString(actualIncrement));
        assertEquals(1.0, actualIncrement[0], 0.000001);
        assertEquals(0.5, actualIncrement[1], 0.000001);
    }

    private GSeries getgSeries(double pointBeingInflunced) throws IOException {
        int fileIdx = findFile(pointBeingInflunced);

        double t0 = gramE12[fileIdx][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = getSavedGSeries(t0, offset);
        return gAtBeta;
    }


}
