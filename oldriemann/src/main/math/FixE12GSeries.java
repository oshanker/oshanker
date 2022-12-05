package math;

import javafx.util.Pair;
import riemann.Interpolate;

import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.LinkedList;

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
    GSeries gAtBeta = null;
    double[] initial = null;
    int midIdxCausingInfluence;
    
    public FixE12GSeries() {
        this(new double[][]{
                {243831.92506103282, -46.745064213360436, -4.265426650034286, 243832.05279114266},
                {243832.1554065314, 81.36514195973275, 17.833911663155945, 243832.3974116106},
                {243832.65617699016, -38.76085656480121, -1.0548384326966262, 243832.71687337887},
                {243832.7750114288, 29.690772236512505, 4.33929464564258, 243832.94901641435},
                {243833.16396249548, -15.05094530616091, -0.5122852558563304, 243833.23656405782},
            },
            //9999
            1999912
        );
    }
    
    public FixE12GSeries(double[][] nextValues, int midIdxCausingInfluence) {
        this.nextValues = nextValues;
        this.midIdxCausingInfluence = midIdxCausingInfluence;
        init();
    }
    
    public void init()  {
        //gAtBeta = getgSeries(pointBeingInflunced[0]);
        try {
            gAtBeta = Interpolate.readGSeries();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public void testIncrementGValuesAtIndices() {
    
        pointBeingInflunced = new double[]{
            nextValues[0][0], nextValues[1][0],
            nextValues[2][0], nextValues[3][0],
            nextValues[4][0],
        };
        initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("initial " + Arrays.toString(initial));

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
    
    public GSeries testChangeToZetaAndDer() {
    
        Pair<double[], double[]> initValues = evaluateNextValues(nextValues, initialPadding, gAtBeta);
        double[] initial = initValues.getKey();
        pointBeingInflunced = initValues.getValue();
        
        System.out.println("initial " + Arrays.toString(initial));
        int[] indices = new int[pointBeingInflunced.length];
        int idxOfCentral = indices.length/2 ;
        for (int i = 0; i < indices.length; i++) {
            indices[i] = midIdxCausingInfluence  - idxOfCentral + i;
        }
        System.out.println(Arrays.toString(indices));
        
        double[][] zetaDerCoeff = changeToZetaAndDer(
            gAtBeta,
            initialPadding,
            pointBeingInflunced,
            initial,
            indices,
            0.125
        );
        double[] neededZetaIncrement = new double[2*pointBeingInflunced.length];
        for (int i = 0; i < pointBeingInflunced.length; i++) {
            int indexIntoNext = i/2;
            if(i%2 == 0) {
                // even i, zero
                if ( pointBeingInflunced[i] != nextValues[indexIntoNext][0] ) {
                    throw new RuntimeException("zero");
                };
                neededZetaIncrement[2 * i] = -initial[2 * i];
                neededZetaIncrement[2 * i + 1] = nextValues[indexIntoNext][1] - initial[2 * i + 1];
            } else {
                // odd i, max
                if ( pointBeingInflunced[i] != nextValues[indexIntoNext][3] ) {
                    System.out.println("pointBeingInflunced[i] != nextValues[indexIntoNext][2]  "
                        + pointBeingInflunced[i] + " != " + nextValues[indexIntoNext][2] );
                    throw new RuntimeException("max");
                };
                neededZetaIncrement[2 * i] = nextValues[indexIntoNext][2]-initial[2 * i];
                neededZetaIncrement[2 * i + 1] =  - initial[2 * i + 1];
            }
        }
    
        LinearEquation linearEquation = new LinearEquation(zetaDerCoeff );
    
        double[] solution = linearEquation.solve(
            neededZetaIncrement
            );
        System.out.println("Required g increment " );
        System.out.println( Arrays.toString(solution));
        gAtBeta.incrementGValuesAtIndices(indices[0], solution);
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        System.out.println("after " );
        System.out.println(Arrays.toString(after));
        double[] actualIncrement = new double[after.length];
        for (int i = 0; i < actualIncrement.length; i++) {
            actualIncrement[i] = after[i] - initial[i];
        }
        System.out.println("actualIncrement " );
        System.out.println( Arrays.toString(actualIncrement));
        
        return gAtBeta;
    }
    
    public static Pair<double[], double[]> evaluateNextValues(
        double[][] nextValues, int initialPadding, GSeries gAtBeta) {
        Pair<double[], double[]> pair;
        LinkedList<Double> values = new LinkedList<>();
        LinkedList<Double> pointBeingInflunced = new LinkedList<>();
        for (int i = 0; i < nextValues.length; i++) {
            pointBeingInflunced.add(nextValues[i][0]);
            values.add(
                gAtBeta.evaluateZeta(nextValues[i][0], initialPadding));
            values.add( gAtBeta.evalDer(
                nextValues[i][0], initialPadding, 0.00025* gAtBeta.spacing));
            if (nextValues[i][3] <= 0) {
                break;
            }
            pointBeingInflunced.add(nextValues[i][3]);
            values.add(
                gAtBeta.evaluateZeta(nextValues[i][3], initialPadding));
            values.add( gAtBeta.evalDer(
                nextValues[i][3], initialPadding, 0.00025* gAtBeta.spacing));
        }
        double[] eval = new double[values.size()];
        for (int i = 0; i < eval.length; i++) {
            eval[i] = values.pollFirst();
        }
        double[] points = new double[pointBeingInflunced.size()];
        for (int i = 0; i < points.length; i++) {
            points[i] = pointBeingInflunced.pollFirst();
        }
        return new Pair<>(eval, points);
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
