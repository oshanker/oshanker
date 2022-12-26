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

import static riemann.StaticMethods.changeToZetaAndDer;
import static riemann.StaticMethods.evaluateAtT;

public class AnalyzeE12GSeries {
    static final int initialPadding = 40;
    private static LinkedList<double[]> zeroInfo;
    private static double maxZeroDev;
    private static double maxDerDev;
    private static double maxMaxDev;
    private static int iMax;
    public static double derepsilon = 1.0E-5;
    double[][] nextValues;
    
    double[] pointBeingInflunced;
    GSeries gAtBeta = null;
    double[] initial = null;
    int midIdxCausingInfluence;
    
    public AnalyzeE12GSeries() {
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
    
    public AnalyzeE12GSeries(
        List<double[]> zeroInfo, int midIdxCausingInfluence, GSeries gAtBeta
    ) {
        this.nextValues = new double[zeroInfo.size()][];
        this.midIdxCausingInfluence = midIdxCausingInfluence;
        this.gAtBeta = gAtBeta;
        for (int row = 0; row < nextValues.length; row++) {
            nextValues[row] = zeroInfo.get(row);
        }
    }
    
    public AnalyzeE12GSeries(double[][] nextValues, int midIdxCausingInfluence) {
        this.nextValues = nextValues;
        this.midIdxCausingInfluence = midIdxCausingInfluence;
        init();
    }
    
    public AnalyzeE12GSeries(double[][] nextValues, int[] midIdxCausingInfluence) {
        this.nextValues = new double[][]{
            {239208.52954627108, 3.4052889402502724, 0.043196050102728534, 239205.96746787973},
            {239208.59376377583, -1.9148244022304233, -0.0942049203052947, 239208.65863886403},
            {239208.72389090285, 2.8047558244728013, 0.23779268540780152, 239208.8225715278},
            {239208.905327034, -5.859201787344891, -0.8499433258858862, 239209.06328053446},
        };
        this.midIdxCausingInfluence = 1961947;
        init();
    }
    
    public void init()  {
        gAtBeta = Interpolate.readGSeries();
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
        double[] actualIncrementInValues = new double[after.length];
        for (int i = 0; i < actualIncrementInValues.length; i++) {
            actualIncrementInValues[i] = after[i] - initial[i];
        }
        System.out.println("actualIncrementInValues " );
        System.out.println( Arrays.toString(actualIncrementInValues));
        
        return gAtBeta;
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
    
    static void resetLimits() {
        maxZeroDev = 3.5;
        maxDerDev = 80.0;
        maxMaxDev = 20.0;
        
    }
    
    public static GSeries testGetSavedGSeries1(
        double firstZero, BufferedReader[] zeroIn, GSeries gAtBeta, int sampleSize
    ) throws IOException {
        iMax = 0;
        int iZero = 0;
        int iDer = 0;
        
        //25 rows zero expectedDer
        
        double[] nextValues = CopyZeroInformation.skipUntil(zeroIn, firstZero);
        if(nextValues[1]<0) nextValues[2] = -nextValues[2];
        int i = 0;
        double z0 = 0, d0 = -1.0, extremumFromFile = -1.0;
    
        zeroInfo = new LinkedList<>();
        //final double stopValue = 243839.5;
        final double stopValue = 243839.0;
        double ZeroDev = 0;
        double DerDev = 0;
        double MaxDev = 0;
    
        for (i = 0; i <= sampleSize; i++) {
            double zeroPosition = nextValues[0];
            double expectedDer =  nextValues[1];
            if (zeroPosition > stopValue) {
                nextValues = CopyZeroInformation.skipUntil(zeroIn, nextValues[0]);
                if (nextValues == null) {
                    System.out.println("End reached " );
                    break;
                }
                if(nextValues[1]<0) nextValues[2] = -nextValues[2];
                System.out.println(i + " End reached " + nextValues[0]);
                break;
            }
            
            double zeta = gAtBeta.evaluateZeta(zeroPosition, initialPadding);
            double der = gAtBeta.evalDer(
                zeroPosition, initialPadding, 0.00025*gAtBeta.spacing);
            double absDer = Math.abs(expectedDer - der);
            boolean maxUpdated = false;
            double absZeta = Math.abs(zeta);
            ZeroDev += absZeta;
            if(absZeta > maxZeroDev){
                maxZeroDev = absZeta;
                iZero = gAtBeta.midIdx;
                maxUpdated = true;
            }
            DerDev += absDer;
            if(absDer > maxDerDev){
                maxDerDev = absDer;
                iDer = gAtBeta.midIdx;
                maxUpdated = true;
            }
            if (maxUpdated) {
                System.out.println("** i " + i + " ===========================");
                System.out.println("zeroPosition " + zeroPosition +
                    " : eval from GSeries: " + zeta + " midIdx " + gAtBeta.midIdx);
                System.out.println(
                    "expectedDer  "
                        + expectedDer
                        + " : eval from GSeries: " + der
                        + " diff(der) " + absDer
                );
            }
            if (i>0) {
                double absMaxDev = findMax(gAtBeta, i, z0, extremumFromFile, zeroPosition);
                MaxDev += absMaxDev;
            }
            double[] zeroEntry = new double[4];
            System.arraycopy(nextValues, 0, zeroEntry, 0, nextValues.length);
            zeroInfo.add(zeroEntry);
            
            z0 = zeroPosition;
            extremumFromFile = nextValues[2];
            if (i < sampleSize) {
                nextValues = CopyZeroInformation.skipUntil(zeroIn, nextValues[0]);
                if(nextValues[1]<0) nextValues[2] = -nextValues[2];
            }
        }
        System.out.println("done");
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
        System.out.println(
            " ZeroDev  " + ZeroDev/i
        );
        System.out.println(
            " DerDev  " + DerDev/i
        );
        System.out.println(
            " MaxDev  " +  MaxDev/i
        );
        
        return gAtBeta;
    }
    
    public static double evalMax(
        GSeries gAtBeta, double x0, double xa, double xb
    ) {
        double positionMax;
        try {
            positionMax = positionMax(gAtBeta, x0, xa, xb);
        } catch (IllegalStateException e) {
            throw e;
        }
        double evalMax = gAtBeta.evaluateZeta(positionMax, initialPadding);
        return evalMax;
    }
    
    public static double positionMax(
        GSeries gAtBeta, double x0, double xa, double xb
    ) {
        double derxb = 0;
        double derxa = gAtBeta.evalDer(
            xa, initialPadding, 0.00025*gAtBeta.spacing);
        if (Math.abs(derxa) < derepsilon) {
            return xa;
        }
        double signumxa = Math.signum(derxa);
        if (signumxa == 0) {
            return xa;
        }
    
        derxb = gAtBeta.evalDer(xb, initialPadding, 0.00025*gAtBeta.spacing);
        if (Math.abs(derxb) < derepsilon) {
            return xb;
        }
        
        double signumxb = Math.signum(derxb);
        if (signumxb == 0) {
            return xb;
        }
        if(signumxa == signumxb) {
            //System.out.println("signumxa == signumxb, " + xa + " " + xb);
            throw new IllegalStateException("signumxa == signumxb");
        }
    
        for (int i = 0; i < 12; i++) {
            if (xb - xa > 0.000001) {
                break;
            }
            double derx0 = gAtBeta.evalDer(x0, initialPadding, 0.00025*gAtBeta.spacing);
            if (Math.abs(derx0) < derepsilon) {
                return x0;
            }
            double signumx0 = Math.signum(derx0);
            if (signumx0 == 0) {
                return x0;
            }
            if (signumx0 == signumxa) {
                xa = x0;
                derxa = gAtBeta.evalDer(
                    xa, initialPadding, 0.00025*gAtBeta.spacing);
                if (Math.abs(derxa) < derepsilon) {
                    return xa;
                }
                signumxa = Math.signum(derxa);
                if (signumxa == 0) {
                    return xa;
                }
            } else {
                xb = x0;
                derxb = gAtBeta.evalDer(xb, initialPadding, 0.00025*gAtBeta.spacing);
                if (Math.abs(derxb) < derepsilon) {
                    return xb;
                }
                signumxb = Math.signum(derxb);
                if (signumxb == 0) {
                    return xb;
                }
            }
            x0 = (xa + xb) / 2;
        }
        return x0;
    }
    
    private static double findMax(
        GSeries gAtBeta, int i, double z0,
        double extremumFromFile,
        double z1
    ) {
        boolean maxMaxUpdated = false;
        double positionMax = positionMax(gAtBeta,(z0+z1)/2, z0, z1);
        double evalMax = gAtBeta.evaluateZeta(positionMax, initialPadding);
        double maxDev = extremumFromFile - evalMax;
        double absMaxDev = Math.abs(maxDev);
        if (absMaxDev > maxMaxDev) {
            maxMaxDev = absMaxDev;
            iMax = gAtBeta.midIdx;
            maxMaxUpdated = true;
        }
        
        if (maxMaxUpdated) {
            System.out.println("** i " + i + " ===========================");
            System.out.println(
                "positionMax " + positionMax
                    + ", eval " + evalMax
                    + " read " + extremumFromFile
                    + " diff(Max) " + maxDev
            );
            System.out.println(" =========*******===========");
        }
        
        double[] oldZero = zeroInfo.getLast();
        oldZero[3] = positionMax;
        return absMaxDev;
    }
    
    private static GSeries getGSeries(double firstZero,  String gbetaSource) throws IOException {
        GSeries gAtBeta;
        System.out.println("** gbetaSource " + gbetaSource + " ======");
        switch (gbetaSource) {
            case "Saved":
                int idx = StaticMethods.findFile(firstZero);
                double t0 = StaticMethods.gramE12[idx][0];
                final BigDecimal offset = BigDecimal.valueOf(1.0E12);
                gAtBeta = StaticMethods.getSavedGSeries(t0, offset);
                break;
            case "Interpolate":
                gAtBeta = Interpolate.readGSeries();
                break;
            default:
                AnalyzeE12GSeries fixE12GSeries = new AnalyzeE12GSeries();
                gAtBeta = fixE12GSeries.testChangeToZetaAndDer();
            
        }
        return gAtBeta;
    }
    
    
    private static void readSavedAndVerifyCross(long first, long last, GSeries gSeries11) {
        double[] cross = new double[]{0, 0};
        int count = 0;
        double base = Interpolate.baseLimit + Interpolate.gramIncr*first;
        double upperLimit = base ;
        long sampleSize = last - first;
        double oldZeta = Double.NEGATIVE_INFINITY;
        while (count < sampleSize  ) {
            double zeta =  Interpolate.evaluateZeta(upperLimit, initialPadding, gSeries11);
            final int nmod2 = count%2;
            if (oldZeta != Double.NEGATIVE_INFINITY) {
                cross[nmod2] += oldZeta * zeta;
            }
            oldZeta = zeta;

            count++;
            upperLimit += Interpolate.gramIncr;
        }
        System.out.println("*** cross Odd " + 2*cross[1]/sampleSize);
        System.out.println("*** cross Even " + 2*cross[0]/sampleSize);
    }
    
    private static void readSavedAndVerify(long first, long last, GSeries gSeries11) {
        int count;
        double[] zetaGramMean = new double[]{0, 0};
        count = 0;
        double base = Interpolate.baseLimit + Interpolate.gramIncr*first;
        double upperLimit = base ;
        long sampleSize = last - first;
        while (count < sampleSize  ) {
            double zeta =  Interpolate.evaluateZeta(upperLimit, initialPadding, gSeries11);
            final int nmod2 = count%2;
            zetaGramMean[nmod2] += zeta;
            count++;
            upperLimit += Interpolate.gramIncr;
        }
        System.out.println("*** zetaGram_MeanOdd " + 2*zetaGramMean[1]/sampleSize);
        System.out.println("*** zetaGram_MeanEven " + 2*zetaGramMean[0]/sampleSize);
    }
    
    public static void prepare(double firstZero, GSeries gAtBeta, int sample) {
        try {
            maxZeroDev = Double.MIN_VALUE;
            maxDerDev = Double.MIN_VALUE;
            maxMaxDev = Double.MIN_VALUE;
            gAtBeta = testGetSavedGSeries1(firstZero, Interpolate.zeroIn, gAtBeta, sample);
            FixE12GSeries fixE12GSeries = new FixE12GSeries(
                zeroInfo.subList(2, 7),
                1999912,
                gAtBeta
            );
        
            gAtBeta = fixE12GSeries.testChangeToZetaAndDer();
            BufferedReader[] zeroIn = null;
            String zerosFile = Rosser.getParam("zerosFile");
            zeroIn = Interpolate.zerosFile(zerosFile);
            maxZeroDev = Double.MIN_VALUE;
            maxDerDev = Double.MIN_VALUE;
            maxMaxDev = Double.MIN_VALUE;
            gAtBeta = testGetSavedGSeries1(firstZero, zeroIn, gAtBeta, sample);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    public static void main(String[] args) {
    
        String gbetaSource = "Interpolate";
        /*
        with poly7, 1.0E-4
maxZeroDev  2.4703315349120514 iZero 1905297
maxDerDev  100.40100690811647 iDer 1985969
maxMaxDev  91.33905303343892 iMax 1961926
with mixed
maxZeroDev  4.523829350361928 iZero 1976284
maxDerDev  120.15046111612506 iDer 1874256
maxMaxDev  33.11496091863871 iMax 1867917
with poly7, 1.0E-6
maxZeroDev  2.4703280886945986 iZero 1905297
maxDerDev  100.4010187254131 iDer 1985969
maxMaxDev  91.34973738048444 iMax 1961926
         */
    
        try {
            GSeries gAtBeta = getGSeries(243831.456494008, gbetaSource);
            //prepare(243831.456494008, gAtBeta, 25);
            //gramVerify(firstZero, gAtBeta);
            //resetLimits();
            maxZeroDev = 2.65;
            maxDerDev = 112.6;
            maxMaxDev = 17;
            /*
246.94434772241235 0.6460803728552249
243840.8746717451 1.5644965586518857

 ZeroDev  0.020492069614571562
 DerDev  0.4955268712260961
 MaxDev  0.08082092773322669
             */
            double firstZero = 248;
            testGetSavedGSeries1(firstZero, Interpolate.zeroIn, gAtBeta,
                2000002);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    private static void gramVerify(double firstZero, GSeries gAtBeta) {
        long firstGramIndex = (long) ((firstZero - Interpolate.baseLimit)/Interpolate.gramIncr);
        long lastGramIndex = (long) ((243839.5956836054 - Interpolate.baseLimit)/Interpolate.gramIncr);
        
        System.out.println(" firstGramIndex " + firstGramIndex );
        System.out.println(" lastGramIndex " + lastGramIndex );
        readSavedAndVerify(firstGramIndex, lastGramIndex, gAtBeta);
        readSavedAndVerifyCross(firstGramIndex, lastGramIndex, gAtBeta);
    }
    
}
