package math;

import javafx.util.Pair;
import riemann.CopyZeroInformation;
import riemann.Interpolate;
import riemann.Rosser;
import riemann.StaticMethods;

import java.io.BufferedReader;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static math.AnalyzeE12GSeries.evalMax;
import static math.AnalyzeE12GSeries.positionMax;
import static riemann.StaticMethods.changeToZetaAndDer;
import static riemann.StaticMethods.evaluateAtT;
import static riemann.StaticMethods.findFile;
import static riemann.StaticMethods.getSavedGSeries;
import static riemann.StaticMethods.gramE12;

public class FixE12GSeries {
    static final int initialPadding = 40;
    public static final double[][] TEST_VALES = {
        {243831.92506103282, -46.745064213360436, -4.265426650034286, 243832.05279114266},
        {243832.1554065314, 81.36514195973275, 17.833911663155945, 243832.3974116106},
        {243832.65617699016, -38.76085656480121, -1.0548384326966262, 243832.71687337887},
        {243832.7750114288, 29.690772236512505, 4.33929464564258, 243832.94901641435},
        {243833.16396249548, -15.05094530616091, -0.5122852558563304, 243833.23656405782},
    };
    
    int midIdxCausingInfluence;
    
    double[] pointBeingInflunced;
    double[][] nextValues;
    GSeries gAtBeta = null;
    
    public FixE12GSeries() {
        this(
            TEST_VALES,
            1999912
        );
    }
    
    public FixE12GSeries(
        List<double[]> zeroInfo, GSeries gAtBeta
    ) {
        this.nextValues = new double[zeroInfo.size()][];
        this.gAtBeta = gAtBeta;
        for (int row = 0; row < nextValues.length; row++) {
            nextValues[row] = zeroInfo.get(row);
        }
    }
    
    public FixE12GSeries(
        List<double[]> zeroInfo, int midIdxCausingInfluence, GSeries gAtBeta
    ) {
        this(zeroInfo, gAtBeta);
        this.midIdxCausingInfluence = midIdxCausingInfluence;
    }
    
    public FixE12GSeries(double[][] nextValues, int midIdxCausingInfluence) {
        this.nextValues = nextValues;
        this.midIdxCausingInfluence = midIdxCausingInfluence;
        init();
    }
    
    public void init()  {
        //gAtBeta = getgSeries(pointBeingInflunced[0]);
        gAtBeta = Interpolate.readGSeries();
    }
    
    
    public double[][] testChangeToZetaAndDer(GSeries gAtBeta) {
    
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
        // this only looks at indices, initial and pointbeinginfluenced.
        // no assumptions regarding zero, max
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
                }
                neededZetaIncrement[2 * i] = -initial[2 * i];
                neededZetaIncrement[2 * i + 1] = nextValues[indexIntoNext][1] - initial[2 * i + 1];
            } else {
                // odd i, max
                if ( pointBeingInflunced[i] != nextValues[indexIntoNext][3] ) {
                    System.out.println("pointBeingInflunced[i] != nextValues[indexIntoNext][2]  "
                        + pointBeingInflunced[i] + " != " + nextValues[indexIntoNext][2] );
                    throw new RuntimeException("max");
                }
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
        double[] actualIncrementInValues = new double[after.length];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            actualIncrementInValues[i] = after[i] - initial[i];
        }
        System.out.println("actualIncrementInValues " );
        System.out.println( Arrays.toString(actualIncrementInValues));
        System.out.println("neededZetaIncrement " );
        System.out.println( Arrays.toString(neededZetaIncrement));
    
        return new double[][]{
            actualIncrementInValues,
            neededZetaIncrement
        };
        
    }
    
    public static Pair<double[], double[]> evaluateNextValues(
        double[][] nextValues, int initialPadding, GSeries gAtBeta) {
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
            // max??
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
    
    public static GSeries testGetSavedGSeries1(
        double firstZero, BufferedReader[] zeroIn, GSeries gAtBeta) throws IOException {
        double maxZeroDev = Double.MIN_VALUE;
        double maxDerDev = Double.MIN_VALUE;
        double maxMaxDev = Double.MIN_VALUE;
        int sampleSize = 25;
    
        int iMax = 0;
        int iZero = 0;
        int iDer = 0;
        
        //25 rows zero expectedDer
        
        double[] nextValues = CopyZeroInformation.skipUntil(zeroIn, firstZero);
        int i = 0;
        double z0 = 0, d0 = -1.0, extremumFromFile = -1.0;
        // cant go below 40
        final int initialPadding = 40;
    
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        
        for (i = 0; i <= sampleSize; i++) {
            double zeroPosition = nextValues[0];
            double expectedDer =  nextValues[1];
            double zeta = gAtBeta.evaluateZeta(zeroPosition, initialPadding);
            double der = gAtBeta.evalDer(
                zeroPosition, initialPadding, 0.00025*gAtBeta.spacing);
            double absDer = Math.abs(expectedDer - der);
            System.out.println("** i " + i + " ===========================");
            boolean maxUpdated = false;
            System.out.println("zeroPosition " + zeroPosition +
                " : eval from GSeries: " + zeta + " midIdx " + gAtBeta.midIdx);
            System.out.println(
                "expectedDer  "
                    + expectedDer
                    + " : eval from GSeries: " + der
                    + " diff(der) " + absDer
            );
            if(Math.abs(zeta) > maxZeroDev){
                maxZeroDev = Math.abs(zeta);
                iZero = i;
                maxUpdated = true;
            }
            if(absDer > maxDerDev){
                maxDerDev = absDer;
                iDer = i;
                maxUpdated = true;
            }
            if (i>0) {
                double positionMax = positionMax(gAtBeta,(z0+zeroPosition)/2, z0, zeroPosition);
                double evalMax = gAtBeta.evaluateZeta(positionMax, initialPadding);
                double maxDev = extremumFromFile - evalMax;
                System.out.println(
                    "positionMax " + positionMax
                        + ", eval " + evalMax
                        + " read " + extremumFromFile
                        + " diff(Max) " + maxDev
                );
                double[] oldZero = zeroInfo.getLast();
                oldZero[3] = positionMax;
                if (oldZero[1] < 0) {
                    oldZero[2] = - oldZero[2];
                }
                if(Math.abs(maxDev) > maxMaxDev){
                    maxMaxDev = Math.abs(maxDev);
                    iMax = i;
                }
                if (maxUpdated) {
                    System.out.println( " =========*******===========");
                } else {
                    System.out.println( " ===========================");
                }
            }
            double[] zeroEntry = new double[4];
            System.arraycopy(nextValues, 0, zeroEntry, 0, nextValues.length);
            zeroInfo.add(zeroEntry);
            System.out.println(
                "nextValues " + Arrays.toString(nextValues)
                    + " extremumFromFile " + extremumFromFile
                    + "; "
            );
            
            z0 = zeroPosition;
            d0 = expectedDer;
            extremumFromFile = d0>0?nextValues[2]:-nextValues[2];
            if(i<sampleSize) {
                nextValues = CopyZeroInformation.skipUntil(zeroIn, nextValues[0]);
            }
        }
        System.out.println("done");
        for (double[] zeroEntry: zeroInfo) {
            System.out.println(
                Arrays.toString(zeroEntry)
            );
        }
        System.out.println(
            "maxZeroDev  " + maxZeroDev
                + " iZero " + iZero
        );
        System.out.println(
            "maxDerDev  " + maxDerDev
                + " iDer " + iDer
        );
        System.out.println(
            "maxMaxDev  " + maxMaxDev
                + " iMax " + iMax
        );
        
        return gAtBeta;
    }
    
    private static GSeries getGSeries(double firstZero,  String gbetaSource) throws IOException {
        GSeries gAtBeta;
        int idx = StaticMethods.findFile(firstZero);
        System.out.println("** gbetaSource " + gbetaSource + " ======");
        switch (gbetaSource) {
            case "Saved":
                double t0 = StaticMethods.gramE12[idx][0];
                final BigDecimal offset = BigDecimal.valueOf(1.0E12);
                gAtBeta = StaticMethods.getSavedGSeries(t0, offset);
                break;
            case "Interpolate":
                gAtBeta = Interpolate.readGSeries();
                break;
            default:
                FixE12GSeries fixE12GSeries = new FixE12GSeries();
                fixE12GSeries.testChangeToZetaAndDer(fixE12GSeries.gAtBeta);
                gAtBeta = fixE12GSeries.gAtBeta;
        }
        return gAtBeta;
    }
    
    public static double[] updateGSeries(
        GSeries gSeries, List<double[]> zeroInfo,
        double[] initial, double[] actual, int[] indices
    ) {
        
        //new double[initialValue.length][2*midIdx_in.length];
        double[][] zetaDerCoeff = gradient(
            gSeries, initial, zeroInfo,
            indices,
            0.125
        );
        double[] neededZetaIncrement = new double[actual.length];
        for (int i = 0; i < actual.length; i++) {
            neededZetaIncrement[i] = actual[i] - initial[i];
        }
        
        LinearEquation linearEquation = new LinearEquation(zetaDerCoeff );
        double determinant = linearEquation.determinant();
        
        double[] solution = linearEquation.solve(
            neededZetaIncrement
        );
        gSeries.incrementGValuesAtIndices(indices[0], solution);
        FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, gSeries);
        double[] initialNoMax =  fixE12GSeries.initialNoMax();
        double[] after = evaluateWithMax(zeroInfo, gSeries, initialNoMax);
        double deviation = getDeviation(actual, after);
        if (deviation > 250) {
            gSeries.decrementGValuesAtIndices(indices[0], solution);
            System.out.println("ignoring =======================");
            LinearEquation.printMatrix(linearEquation.coefficients);
            System.out.println("indices " + Arrays.toString(indices));
            System.out.println("after " + Arrays.toString(after));
            System.out.println("          =======================");
            return new double[] {getDeviation(actual, initial), determinant};
        }
    
        return new double[] {deviation, determinant};
    }
    
    private static double getDeviation(double[] actual, double[] after) {
        double deviation = 0;
        int sample = 0;
        for (int i = 0; i < after.length; i++) {
            if (i%3 == 1) {
                //ignore derivative
                continue;
            }
            sample++;
            deviation += Math.abs(after[i] - actual[i]);
        }
        deviation /= sample;
        return deviation;
    }
    
    public static boolean advanceZeroInfo(
        BufferedReader[] zeroIn, double firstZero, LinkedList<double[]> zeroInfo
    ) {
        if (zeroInfo.get(1)[0] < firstZero) {
            double[] nextValues = CopyZeroInformation.skipUntil(zeroIn, 0);
            if(nextValues == null) {
                return false;
            }
            if(nextValues[1]<0){
                nextValues[2]=-nextValues[2];
            }
            zeroInfo.pollFirst();
            zeroInfo.add(nextValues);
            return true;
        }
        return false;
    }
    
    public static void initZeroInfo(
        BufferedReader[] zeroIn, double firstZero,
        LinkedList<double[]> zeroInfo, int desiredSize
    ){
        if (zeroInfo.size() < desiredSize) {
            int needed = desiredSize - zeroInfo.size();
            for (int i = 0; i < needed; i++) {
                double[] nextValues = CopyZeroInformation.skipUntil(zeroIn, firstZero);
                if(nextValues[1]<0){
                    nextValues[2]=-nextValues[2];
                }
                zeroInfo.add(nextValues);
            }
        }
    }
    
    static void printZeroInfo(LinkedList<double[]> zeroInfo) {
        for (double[] doubles : zeroInfo) {
            System.out.println(Arrays.toString(doubles));
        }
    }
    
    static double[][] printZeroInfoWithMax(GSeries gAtBeta, LinkedList<double[]> zeroInfo) {
        int size =  zeroInfo.size();
        double[] initial = new double[3*size-1];
        double[] actual = new double[3*size-1];
        // initial, actual
        int midIdxCausingInfluence = -1;
        double[] oldzero = null;
        for (int i = 0; i < size; i++) {
            double[] zero = zeroInfo.get(i);
            initial[3*i] = gAtBeta.evaluateZeta(zero[0], initialPadding);
            if(midIdxCausingInfluence < 0){
                midIdxCausingInfluence = gAtBeta.midIdx;
                System.out.println("*** gAtBeta.midIdx " + midIdxCausingInfluence
                    + " " + (zero[0]));
        
            }
            initial[3*i+1] = gAtBeta.evalDer(
                zero[0], initialPadding, 0.00025* gAtBeta.spacing);
            actual[3*i] = 0;
            actual[3*i+1] = zero[1];
            if (i>0) {
                double positionMax = positionMax(gAtBeta,
                    (oldzero[0] + zero[0]) / 2, oldzero[0], zero[0]);
                double evalMax = gAtBeta.evaluateZeta(positionMax, initialPadding);
                initial[3 * i - 1] = evalMax;
                actual[3 * i - 1] = oldzero[2];
                System.out.println(positionMax + " " + evalMax);
            }
            oldzero = zero;
            System.out.println(Arrays.toString(zero));
        }
        
        // indices
        int[] indices = new int[(3*size-1)/2];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = midIdxCausingInfluence + i;
        }
        
        // update
        double[] deviation = updateGSeries(
            gAtBeta, zeroInfo,
            initial,  actual,  indices
        );
        System.out.println("deviation " + deviation[0]
            + " det " + deviation[1]);
        double[][] ret = new double[][] {
            deviation, new double[]{midIdxCausingInfluence}};
        return ret;
    }
    
    static double[] applyFix(
        GSeries gAtBeta, int midIdxCausingInfluence, LinkedList<double[]> zeroInfo
    ) {
        int size = zeroInfo.size();
        FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, gAtBeta);
    
        double[] initialNoMax =  fixE12GSeries.initialNoMax();
        boolean useNoMax = false;
        for (int i = 1; i < initialNoMax.length-1; i+=2) {
            double signumxa = Math.signum(initialNoMax[i]);
            double signumxb = Math.signum(initialNoMax[i+2]);
            if (
                signumxa == 0 ||
                signumxb == 0 ||
                signumxa == signumxb
            ) {
                useNoMax = true;
                break;
            }
        }
        if(useNoMax) {
            double[][] ret = fixE12GSeries.testChangeToZetaAndDerNoMax(
                fixE12GSeries.gAtBeta, midIdxCausingInfluence, false);
            double[] deviation = ret[2];
            return deviation;
        }
        double[] initial = evaluateWithMax(zeroInfo, gAtBeta, initialNoMax);
        double[] actual = actual(zeroInfo);
        
        // indices
        int[] indices = new int[(3*size-1)/2];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = midIdxCausingInfluence + i;
        }
//        double deviation = getDeviation(actual, initial);
//        if (deviation < 0.001) {
//            return new double[] {deviation, 1000};
//        }
        
        // update
        double[] deviationDet = updateGSeries(
            gAtBeta, zeroInfo,
            initial,  actual,  indices
        );
        if (deviationDet[0] > 250) {
                System.out.println("actual " + Arrays.toString(actual));
                System.out.println("initial " + Arrays.toString(initial));
        }
        return deviationDet;
    }
    
    public static double[][] gradient(
        GSeries gSeries, double[] initialValue, List<double[]> zeroInfo,
                int[] midIdx_in,
               double increment
    ) {
        int size = zeroInfo.size();
        double[] oldzero = null;
        
        double[][] gradient = new double[initialValue.length][2*midIdx_in.length];
        for (int gSeriesIndex = 0; gSeriesIndex < midIdx_in.length; gSeriesIndex++) {
            int midIdx = midIdx_in[gSeriesIndex];
            gSeries.incrementGValueAtIndex(midIdx, new double[]{increment, 0});
            oldzero = null;
            for (int idxZeroInfo = 0; idxZeroInfo < size; idxZeroInfo++) {
                double[] zero = zeroInfo.get(idxZeroInfo);
                gradient[3*idxZeroInfo][2*gSeriesIndex] = (gSeries.evaluateZeta(zero[0], initialPadding)-initialValue[3*idxZeroInfo])/ increment;
                double a0ValAtZero = gSeries.evalDer(
                    zero[0], initialPadding, 0.00025 * gSeries.spacing);
                gradient[3*idxZeroInfo+1][2*gSeriesIndex] = (a0ValAtZero - initialValue[3*idxZeroInfo + 1]) / increment;
                if (idxZeroInfo > 0) {
                    gradient[3 * idxZeroInfo - 1][2 * gSeriesIndex] = (evalMax(gSeries,
                        (oldzero[0] + zero[0]) / 2, oldzero[0], zero[0]) - initialValue[3*idxZeroInfo - 1]) / increment;
                }
                oldzero = zero;
            }
            
            gSeries.incrementGValueAtIndex(midIdx, new double[]{-increment, 0});
            
            // change second index
            
            gSeries.incrementGValueAtIndex(midIdx, new double[]{0, increment});
    
            oldzero = null;
            for (int idxZeroInfo = 0; idxZeroInfo < size; idxZeroInfo++) {
                double[] zero = zeroInfo.get(idxZeroInfo);
                gradient[3*idxZeroInfo][2*gSeriesIndex+1] = (gSeries.evaluateZeta(zero[0], initialPadding)-initialValue[3*idxZeroInfo])/ increment;
    
                double a0ValAtZero = gSeries.evalDer(
                    zero[0], initialPadding, 0.00025 * gSeries.spacing);
                gradient[3*idxZeroInfo+1][2*gSeriesIndex+1] = (a0ValAtZero - initialValue[3*idxZeroInfo + 1]) / increment;
                if (idxZeroInfo > 0) {
                    gradient[3 * idxZeroInfo - 1][2 * gSeriesIndex + 1] = (evalMax(gSeries,
                        (oldzero[0] + zero[0]) / 2, oldzero[0], zero[0]) - initialValue[3*idxZeroInfo - 1]) / increment;
                }
                oldzero = zero;
            }
            gSeries.incrementGValueAtIndex(midIdx, new double[]{0, -increment});
        }
        
        return gradient;
    }
    
    public static double[] actual(
        List<double[]> zeroInfo
    ) {
        int size = zeroInfo.size();
        double[] ret = new double[3* size -1];
        for (int i = 0; i < size; i++) {
            double[] zero = zeroInfo.get(i);
            ret[3*i] = 0;
            ret[3*i+1] = zero[1];
            if (i < size-1) {
                ret[3 * i + 2] = zero[2];
            }
        }
        return ret;
    }
   
    public static double[] evaluateWithMax(
        List<double[]> zeroInfo,  GSeries gAtBeta, double[] initialNoMax
    ) {
        double[] ret = new double[3*zeroInfo.size()-1];
        double[] oldzero = null;
        for (int i = 0; i < zeroInfo.size(); i++) {
            double[] zero = zeroInfo.get(i);
            if (i>0) {
                double evalMax = evalMax(gAtBeta,
                    (oldzero[0] + zero[0]) / 2, oldzero[0], zero[0]);
                ret[3 * i - 1] = evalMax;
            }
            oldzero = zero;
            ret[3*i] = initialNoMax[2*i];
            ret[3*i+1] = initialNoMax[2*i+1];
        }
        return ret;
    }
    
    public static void main(String[] args) throws IOException {
        fixGSeries01();
    }
    
    double[] initialNoMax(){
        double[] initial =  new double[2*nextValues.length];
        for (int i = 0; i < nextValues.length; i++) {
            initial[2*i] =
                gAtBeta.evaluateZeta(nextValues[i][0], initialPadding);
            initial[2*i+1] =
                gAtBeta.evalDer(
                    nextValues[i][0], initialPadding, 0.00025 * gAtBeta.spacing);
        }
        return initial;
    }
    
    public double[][] testChangeToZetaAndDerNoMax(
        GSeries gAtBeta, int idxCausingInfluence, boolean verbose) {
        
        pointBeingInflunced = new double[nextValues.length];
        int[] indices = new int[pointBeingInflunced.length];
        double[] initial =  initialNoMax();
        
        for (int i = 0; i < nextValues.length; i++) {
            pointBeingInflunced[i] = nextValues[i][0];
            if (idxCausingInfluence < 0) {
                idxCausingInfluence = gAtBeta.midIdx;
            }
            indices[i] = idxCausingInfluence + i;
        }
        if (verbose) {
            System.out.println("initial " + Arrays.toString(initial));
            System.out.println(Arrays.toString(indices));
        }
        // this only looks at indices, initial and pointbeinginfluenced.
        // no assumptions regarding zero, max
        double[][] zetaDerCoeff = changeToZetaAndDer(
            gAtBeta,
            initialPadding,
            pointBeingInflunced,
            initial,
            indices,
            0.125
        );
        double[] actual = new double[2*pointBeingInflunced.length];
        double[] neededZetaIncrement = new double[2*pointBeingInflunced.length];
        for (int i = 0; i < pointBeingInflunced.length; i++) {
            neededZetaIncrement[2 * i] = -initial[2 * i];
            neededZetaIncrement[2 * i + 1] = nextValues[i][1] - initial[2 * i + 1];
            actual[2 * i] = 0;
            actual[2 * i + 1] = nextValues[i][1];
        }
        
        LinearEquation linearEquation = new LinearEquation(zetaDerCoeff );
        double determinant = linearEquation.determinant();
        
        double[] solution = linearEquation.solve(
            neededZetaIncrement
        );
        if (verbose) {
            System.out.println("Required g increment ");
            System.out.println(Arrays.toString(solution));
        }
        
        gAtBeta.incrementGValuesAtIndices(indices[0], solution);
        
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        if (verbose) {
            System.out.println("after ");
            System.out.println(Arrays.toString(after));
        }
        double[] actualIncrementInValues = new double[after.length];
        double deviation = 0;
        int sample = 0;
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            actualIncrementInValues[i] = after[i] - initial[i];
            if (i%2 == 0) {
                sample++;
                deviation += Math.abs(after[i] - actual[i]);
            }
        }
        deviation /= sample;
        if (verbose) {
            System.out.println("actualIncrementInValues ");
            System.out.println(Arrays.toString(actualIncrementInValues));
            System.out.println("neededZetaIncrement ");
            System.out.println(Arrays.toString(neededZetaIncrement));
        }
        
        return new double[][]{
            actualIncrementInValues,
            neededZetaIncrement,
            new double[] {deviation, determinant}
        };
    }
    
    private static void fixGSeries01() {
        GSeries gAtBeta = Interpolate.readGSeries();
        int R = gAtBeta.gAtBeta.length;
        double begin = gAtBeta.begin + (initialPadding-18)*gAtBeta.spacing;
        System.out.println(begin + " " + gAtBeta.evaluateZeta(begin+ gAtBeta.spacing/2, initialPadding));
        System.out.println("gAtBeta.midIdx " + gAtBeta.midIdx
         + " " + (gAtBeta.begin + gAtBeta.midIdx*gAtBeta.spacing));
        System.out.println("==========");
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        int desiredSize = 3;
        initZeroInfo(
            Interpolate.zeroIn, begin + 6*gAtBeta.spacing, zeroInfo, desiredSize);
        System.out.println("==========");
        double[][] initialRet = printZeroInfoWithMax(gAtBeta, zeroInfo);
        int midIdxCausingInfluence = (int)initialRet[1][0];
        System.out.println("==========");
        int devCount = 0;
        int worseCount = 0;
        double maxDev = Double.MIN_VALUE;
        double minDet = Double.MAX_VALUE;
        for (int iter = 0; iter < 50000; iter++) {
            midIdxCausingInfluence++;
            double nextValue = gAtBeta.begin + midIdxCausingInfluence*gAtBeta.spacing;
            advanceZeroInfo(Interpolate.zeroIn, nextValue, zeroInfo);
            double[] ret = null;
            boolean useMax = true;
            try {
                ret = applyFix(gAtBeta, midIdxCausingInfluence, zeroInfo);
            } catch (IllegalStateException e) {
                // do a fit with only zeros and derivative
                //System.out.println(" skipping, xa xb same sign, " + midIdxCausingInfluence);
                useMax = false;
                continue;
//                FixE12GSeries fixE12GSeries = new FixE12GSeries(zeroInfo, gAtBeta);
//                double[][] ret1 = fixE12GSeries.testChangeToZetaAndDerNoMax(gAtBeta, midIdxCausingInfluence);
//                ret = ret1[2];
            }
            
            double deviation = ret[0];
            double det = Math.abs(ret[1]);
            
//            System.out.println(midIdxCausingInfluence + " **** " + deviation
//                + ", det " + det);
            
            if (deviation > 0.01) {
                if (deviation > 1.0) {
                    System.out.println("midIdx " + midIdxCausingInfluence +
                        " nextValue " + nextValue);
                    System.out.println("deviation too large, " + deviation
                        + ", det " + det);
                    System.out.println("=====** " + ++devCount);
                    if (deviation > 250) {
                        for (int zeroIdx = 0; zeroIdx < zeroInfo.size(); zeroIdx++) {
                            System.out.println(" " + Arrays.toString(zeroInfo.get(zeroIdx)));
                        }
                        throw new IllegalStateException("check this value!");
                    }
                }
                if (useMax) {
                    double oldDeviation = deviation;
                    ret = applyFix(gAtBeta, midIdxCausingInfluence, zeroInfo);
                    // need signum check here too.
                    deviation = ret[0];
                    det = Math.abs(ret[1]);
                    if (deviation > 0.1) {
                        System.out.println("midIdx " + midIdxCausingInfluence +
                            " nextValue " + nextValue + " oldDeviation " + oldDeviation);
                        System.out.println("deviation too large (after retry), " + deviation);
                        if (deviation > oldDeviation) {
                            System.out.println("=====!!!!!!!!!=== " + ++worseCount);
                        }
                        System.out.println("===============");
                    }
                }
            }
            if (deviation > maxDev) {
                maxDev = deviation;
            }
            if (det == 0) {
                System.out.println(" det == 0, " + midIdxCausingInfluence);
            }
            if (det < minDet) {
                minDet = det;
            }
        }
        System.out.println("maxDev " + maxDev + " minDet " + minDet);
    
        System.out.println("==========");
        double nextValue = gAtBeta.begin + midIdxCausingInfluence*gAtBeta.spacing;
        System.out.println("Final: midIdx " + midIdxCausingInfluence + " nextValue " + nextValue);
        System.out.println("==========");
        double end = gAtBeta.begin + (R-22)*gAtBeta.spacing;
        System.out.println(end + " " + gAtBeta.evaluateZeta(end - gAtBeta.spacing/2, initialPadding));
        System.out.println("gAtBeta.midIdx " + gAtBeta.midIdx
            + " " + (gAtBeta.begin + gAtBeta.midIdx*gAtBeta.spacing));
    }
    
    private static void oldMain() {
        double firstZero = 243831.456494008;
        String gbetaSource = "Interpolate";
        
        try {
            GSeries gAtBeta = getGSeries(firstZero, gbetaSource);
            gAtBeta = testGetSavedGSeries1(firstZero, Interpolate.zeroIn, gAtBeta);
            LinkedList<double[]> zeroInfo = new LinkedList<>();
            FixE12GSeries fixE12GSeries = new FixE12GSeries(
                zeroInfo.subList(2, 7),
                1999912,
                gAtBeta
            );
            
            fixE12GSeries.testChangeToZetaAndDer(fixE12GSeries.gAtBeta);
            gAtBeta = fixE12GSeries.gAtBeta;
            BufferedReader[] zeroIn = null;
            String zerosFile = Rosser.getParam("zerosFile");
            zeroIn = Interpolate.zerosFile(zerosFile);
            gAtBeta = testGetSavedGSeries1(firstZero, zeroIn, gAtBeta);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
}
