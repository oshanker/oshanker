package math;

import riemann.Interpolate;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;

import static riemann.StaticMethods.changeToDer;
import static riemann.StaticMethods.changeToZeta;
import static riemann.StaticMethods.changeToZetaAndDer;
import static riemann.StaticMethods.evaluateAtT;
import static riemann.StaticMethods.findFile;
import static riemann.StaticMethods.getSavedGSeries;
import static riemann.StaticMethods.gramE12;

public class FixE12GSeries {
    static final int initialPadding = 40;
    double[][] nextValues;
    
    double[] pointBeingInflunced;
    static GSeries gAtBeta = null;
    static double[] initial = null;
    int midIdxCausingInfluence;
    
    public FixE12GSeries() {
        this(new double[][]{
                {243831.92506103282, -46.745064213360436, 4.265426650034286},
                {243832.1554065314, 81.36514195973275, 17.833911663155945},
                {243832.65617699016, -38.76085656480121, 1.0548384326966262},
                {243832.7750114288, 29.690772236512505, 4.33929464564258},
                {243833.16396249548, -15.05094530616091, 0.5122852558563304},
            },
            //9999
            1999912
        );
    }
    
    public FixE12GSeries(double[][] nextValues, int midIdxCausingInfluence) {
        this.nextValues = nextValues;
        this.midIdxCausingInfluence = midIdxCausingInfluence;
        pointBeingInflunced = new double[]{
            nextValues[0][0], nextValues[1][0],
            nextValues[2][0], nextValues[3][0],
            nextValues[4][0],
        };
        init();
    }
    
    public void init()  {
        //gAtBeta = getgSeries(pointBeingInflunced[0]);
        try {
            gAtBeta = Interpolate.readGSeries();
        } catch (IOException e) {
            e.printStackTrace();
        }
    
        initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("initial " + Arrays.toString(initial));
    }
    
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
    }
    
    public void testChangeToZetaAndDer() {
        int[] indices = {
            midIdxCausingInfluence-2, midIdxCausingInfluence-1, midIdxCausingInfluence,
            midIdxCausingInfluence+1, midIdxCausingInfluence+2,
        };
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
            new double[]{
                //0.029814082487805593, 0.1994990926788418, 0.08136282505069659, 2.125999647229591, 0.9999999999999996, 0.5000000000951985, 0.38281568257124876, -7.700267711652717,
                -1.5533670572054348E-7, 3.455601111568285E-5, 0.002939713190428961, 0.17302076126421184, 1.0000000000000138,
                0.5000000000195968, -0.6319767951186775, 5.067432651792831,
                -4.3272877404021415E-5, 0.0026806606801486055

            });
        System.out.println("Required g increment " );
        System.out.println( Arrays.toString(solution));
        gAtBeta.incrementGValuesAtIndices(indices[0], solution);
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("after " + Arrays.toString(after));
        double[] actualIncrement = new double[after.length];
        for (int i = 0; i < actualIncrement.length; i++) {
            actualIncrement[i] = after[i] - initial[i];
        }
        System.out.println("actualIncrement " );
        System.out.println( Arrays.toString(actualIncrement));


//        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
//        System.out.println("after " );
//        System.out.println( Arrays.toString(after));
    }
    
    private static GSeries getgSeries(double pointBeingInflunced)  {
        int fileIdx = findFile(pointBeingInflunced);

        double t0 = gramE12[fileIdx][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = null;
        try {
            gAtBeta = getSavedGSeries(t0, offset);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return gAtBeta;
    }

    public static void main(String[] args) {
        FixE12GSeries fixE12GSeries = new FixE12GSeries();
        fixE12GSeries.testChangeToZetaAndDer();
        //fixE12GSeries.testIncrementGValuesAtIndices();
    }

}
