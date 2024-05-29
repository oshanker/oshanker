/**
 *
 */
package math;

import static math.AnalyzeE12GSeries.positionMax;
import static math.FixGSeries.evaluateAtT;
import static org.junit.Assert.*;
import static riemann.StaticMethods.findFile;
import static riemann.StaticMethods.getSavedGSeries;
import static riemann.StaticMethods.gramE12;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.LinkedList;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import riemann.*;

/**
 * @author oshanker
 *
 */
public class GSeriesTest {
    static final int zRange = 39;
    final double epsilon = 0.0005;
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(8);
        nf.setMaximumFractionDigits(8);
        nf.setGroupingUsed(false);
    }

    /**
     * @throws java.lang.Exception
     */
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        Gram.initLogVals((int)(398942/2.2));
    }

    /**
     * @throws java.lang.Exception
     */
    @AfterClass
    public static void tearDownAfterClass() throws Exception {
    }

    /**
     * @throws java.lang.Exception
     */
    @Before
    public void setUp() throws Exception {
    }

    /**
     * @throws java.lang.Exception
     */
    @After
    public void tearDown() throws Exception {
    }

    @Test
    public void testIncrementGValuesAtIndices() throws IOException {
        double[] nextValues = new double[]{
                243831.456494008, -22.69554476177354, 1.538114456203189};
        double pointBeingInflunced = nextValues[0];
        int fileIdx = findFile(pointBeingInflunced);
        final int initialPadding = 40;

        double t0 = gramE12[fileIdx][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = getSavedGSeries(t0, offset);
        double[] initial = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        int midIdxCausingInfluence  = gAtBeta.midIdx;
        double[][] requiredGIncrements = {
                {0.13861482493228286, 0.8115131385779506},
                {0, 0}
        };
        gAtBeta.incrementGValuesAtIndices(midIdxCausingInfluence, requiredGIncrements);
        double[] after = evaluateAtT(pointBeingInflunced, initialPadding, gAtBeta);
        double[] actualIncrement = new double[after.length];
        for (int i = 0; i < actualIncrement.length; i++) {
            actualIncrement[i] = after[i] - initial[i];
        }
        System.out.println("actualIncrement " + Arrays.toString(actualIncrement));
        assertEquals(1.0, actualIncrement[0], 0.000001);
        assertEquals(0.5, actualIncrement[1], 0.000001);
    }

    /**
     * Test method for t = 2.7E9.
     */
    @Test @Ignore
    public void testE27E9() {
        int k0 = 1, k1=206393;
        int R = 10000;
        long init= System.currentTimeMillis();
        BigDecimal offset = BigDecimal.valueOf(267653395647L);
        double begin = 1.87383225;
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
        GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, R);
        long end = System.currentTimeMillis();
        System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
        System.out.println(gAtBeta.riemannZeta(gAtBeta.gAtBeta[0], begin));
        //g  : [-0.33143958775035764, 0.0733285174786178] 1287.5146091794
        double zero = 1287.5146091794;
        double[] gFromBLFI = gAtBeta.blfiSumWithOffsetSmallT( zero, 4);
        double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
        System.out.println("g  : " + Arrays.toString(gFromBLFI) + " zeta " + zeta);
        assertTrue(Math.abs(gFromBLFI[0] - (-0.33143958775035764)) + Math.abs(gFromBLFI[1] - 0.0733285174786178) < 0.000005);
    }

    @Test @Ignore
    public void testGenerateGramE12() throws Exception{
        int R = 30040;
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        double begin = Gram.gram(offset, 243.77756012466054 );
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin + (R-40)*0.2436/2)/(2*Math.PI)));
        String zerosFile = "/Users/shankero/Documents/tmp/1e12.zeros.1001_10001002";
        BufferedReader[] zeroIn = {new BufferedReader(new FileReader(zerosFile ))};
        double upperLimit = begin + (R-40)*incr;
        Rosser.ZeroInfo zeroInput = Rosser.readZeros(upperLimit, null, zeroIn,  null);
        System.out.println("{" +begin + ", " + zeroInput.lastZero[0]+ "},");
        for (int i = 0; i < 34; i++) {
            begin = upperLimit;
            incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin + (R-40)*incr/2)/(2*Math.PI)));
            upperLimit = begin + (R-40)*incr;
            zeroInput = Rosser.readZeros(upperLimit, null, zeroIn,  zeroInput.nextValues);
            System.out.println( "{" + begin + ", " + zeroInput.lastZero[0] + "},");
        }
    }

    /**
     * Test method for calculating calculating and storing G at E12.
     */
    @Test  @Ignore
    public void testStoreGE12() throws Exception{
        for (int i = 0; i < StaticMethods.gramE12.length; i++) {
            storeGE12(StaticMethods.gramE12[i][1], StaticMethods.gramE12[i][0]);
        }
    }

    @Test @Ignore
    public void testCorrelationE12() throws Exception{
        nf.setMinimumFractionDigits(3);
        nf.setMaximumFractionDigits(3);
        nf.setGroupingUsed(false);
        File file = new File("out/gzetaCorrelation/gzeta8.csv");
        if (!file.getParentFile().exists()) {
            file.getParentFile().mkdirs();
        }
        PrintWriter out = new PrintWriter(file);
        //PrintWriter out = new PrintWriter(System.out);
        final int k = 8;
        final double[][] prod = new double[2*k][2*k];
        int size = 0;
        for (int i = 0; i < StaticMethods.gramE12.length; i++) {
            size += testCorrelationE12(i, out, prod );
        }
        for (int i = 0; i < 2*k; i++) {
            for (int j = 0; j < 2*k; j++) {
                if(j>0){out.print(", ");}
                out.print(nf.format(prod[i][j]/size));
            }
            out.println();
        }
        out.close();
    }

    @Test @Ignore
    // need  out/gzetaE12/calcHist12.csv (No such file or directory)
    public void testSplineFit() throws Exception{
        //needs fixing
        int k = 12;
        File file = new File("out/gzetaE12/calcHist" + k +  ".csv");
        BufferedReader zeroIn = new BufferedReader(new FileReader(file));
        int k2 = 2*k;
        double[][] vals = new double[zRange-2][k2];
        int i = 0;
        double[] x = new double[zRange-2];
        double[] y= new double[zRange-2];
        double[] y4= new double[zRange-2];
        double[] y2= new double[zRange-2];
        String in = zeroIn.readLine();
        for (i = 0; i < 2; i++) {
            in = zeroIn.readLine();
        }
        int iZero = zRange/2-1;
        int phiIndex = 2;
        double phi = phiIndex*Math.PI/12.0;
        for (i = 0; i < zRange-2; i++) {
            String[] line = in.split(",");
            x[i]= Double.parseDouble(line[0].trim());
            int j = phiIndex + 1;
            y[i]= Double.parseDouble(line[j].trim());
            y4[i]= Double.parseDouble(line[4].trim());
            y2[i]= Double.parseDouble(line[7].trim());
            for (j = 1; j <= vals[i].length; j++) {
                vals[i][j-1]= Double.parseDouble(line[j].trim());
            }

            in = zeroIn.readLine();
        }
        assertEquals(0, x[iZero], 0.000001);
        zeroIn.close();
        NormalizedSpline normalizedSpline = new NormalizedSpline(x[0], x[x.length-1], y);
        NormalizedSpline spline4 = new NormalizedSpline(x[0], x[x.length-1], y4);
        NormalizedSpline spline2 = new NormalizedSpline(x[0], x[x.length-1], y2);
        double mult = vals[iZero][phiIndex]/(y4[iZero]*y2[iZero]);
        assertEquals(0.33099349, y4[iZero], 0.000001);
        assertEquals(0.34465536, y2[iZero], 0.000001);
        assertEquals(0.33099349, vals[iZero][3], 0.000001);
        assertEquals(0.34465536, vals[iZero][6], 0.000001);
        System.out.println(phi*180/Math.PI + " check " + mult*0.33099349*0.34465536);
        final double root2 = Math.sqrt(2.0);

        final double cosphi = Math.cos(phi);
        final double cos2phi = Math.cos(2*phi);
        for(i = 10; i<26;i++) {
            double xi = x[i];
//			double pred = mult*spline4.eval(xi*root2*cosphi)
//        			*spline2.eval(xi*cos2phi);
            double pred = (1+cos2phi)*spline4.eval(xi*root2*cosphi)
                    -cos2phi*spline2.eval(xi*cos2phi);

            System.out.println(nf.format(xi) + " "
                    + nf.format(vals[i][phiIndex]) + " "
                    // + nf.format(spline2.eval(xi*Math.cos(2*phi))) + " "
                    + nf.format(pred) + " "
                    + nf.format(Math.abs(vals[i][phiIndex]-pred)) + " "
            );
        }
    }

    @Test //@Ignore  
    public void testSymmetryRelations() throws Exception{
        //check the symmetry and antisymmetry relations from the output of distributions
        int k = 12;
        //testWriteZetaPhiE12
        File file = new File("out/gzetaE12/calcHist_fine" + k +  ".csv");
        BufferedReader zeroIn = new BufferedReader(new FileReader(file));
        int k2 = 2*k;
        double[][] vals = new double[zRange][k2];
        String in = zeroIn.readLine();
        in = zeroIn.readLine();
        int i = 0;
        while (in != null) {
            String[] line = in.split(",");
            for (int j = 1; j <= vals[i].length; j++) {
                vals[i][j-1]= Double.parseDouble(line[j].trim());
            }
            in = zeroIn.readLine();
            i++;
        }
        zeroIn.close();
        //antisymmetry
        double maxdiff = 0;
        for (i = 0; i < zRange; i++) {
            for (int j = 0; j < k; j++) {
                double diff = -100;
                final double other = vals[zRange-1-i][k+j];
                diff = Math.abs( vals[i][j] - other);
                if(diff> maxdiff) {maxdiff =diff;}
                System.out.print( j + " " + nf.format(vals[i][j]) + " "
                        + nf.format(other) + " " + nf.format(diff) + " | ");
                //System.out.print(i + " " + j + " " + nf.format(diff) + " ");
                //assertEquals("", 0, diff, 0.003);
                //temp
                assertEquals("", 0, diff, 0.03);
            }
            System.out.println(" | " );
        }
        System.out.println(" | " + maxdiff);
        //symmetry
        maxdiff = 0;
        for (i = 0; i < zRange; i++) {
            for (int j = 1; j < k; j++) {
                double diff = -100;
                final double other = vals[i][k2-j];
                diff = Math.abs( vals[i][j] - other);
                if(diff> maxdiff) {maxdiff =diff;}
                System.out.print( j + " " + nf.format(vals[i][j]) + " "
                        + nf.format(other) + " " + nf.format(diff) + " | ");
                //System.out.print(i + " " + j + " " + nf.format(diff) + " ");
                //assertEquals("", 0, diff, 0.0035);
                //temp
                assertEquals("", 0, diff, 0.03);
            }
            System.out.println(" | " );
        }
        System.out.println(" | " + maxdiff);
    }

    @SuppressWarnings("unused")
    @Test //@Ignore
    public void testWriteZetaPhiE12() throws Exception{
        //zetaQuantile.R

        //calculate Z, mean(Z), store Z
        //z, quantile, and hist
        //This could provide input for zetaHist.R
        // after that is run, use the testSymmetryRelations method (MoreGseriesTest)
        //zetaQuantile.R or oshanker/python/riemann/plot_distribution.py

        final int k = 12;
        final File baseDir = new File("out/gzetaE12/");
        final File outputZFile = new File(baseDir,
                "gzeta_calc" + k + ".csv");
        //      "gzeta_calc" + k + ".csv");
        if (!outputZFile.getParentFile().exists()) {
            outputZFile.getParentFile().mkdirs();
        }
        //final PrintWriter outputZ = new PrintWriter(outputZFile);
        final PrintWriter outputZ = null;

        final double[][] gramSum = new double[1][2*k];
        final Histogram[] hist = new Histogram[2*k];
//        final double min = -9.25;
//        final double max = 9.25;
        final double min = -.925;
        final double max = 0.925;
        final int binCount = zRange-2;
        for (int i = 0; i < hist.length; i++) {
            hist[i] = new Histogram(min, max, binCount);
        }

        //if we want cross-product
        //final double[][] gramSum = new double[2][2*k];

        for (int i = 0; i < StaticMethods.gramE12.length; i++) {
            double[][] cross = new double[gramSum.length][2*k];
            writeZetaPhi(i, outputZ, cross, hist);
            for (int j = 0; j < 2*k; j++) {
                gramSum[0][j] += cross[0][j];
                for (int i1 = 1; i1 < gramSum.length; i1++) {
                    gramSum[i1][j] += cross[i1][j];
                }
            }
//	        System.out.println(nf.format(gramSum[1][0]/(3.154*(i+1))) +
//	        		"\t " + nf.format(gramSum[1][k]/(3.154*(i+1))));
        }

        if(outputZ != null) {
            outputZ.close();
        }

        final double crossNorm = 2*(1.577)* StaticMethods.gramE12.length;
        for (int j = 0; j < 2*k; j++) {
            gramSum[0][j] /= StaticMethods.gramE12.length;
            assertEquals(2.00*Math.cos(j*Math.PI/k), gramSum[0][j], 0.001);
            for (int i = 1; i < gramSum.length; i++) {
                gramSum[i][j] /= crossNorm;
            }
        }
        for (int i = 0; i < gramSum.length; i++) {
            for (int j = 0; j < 2*k; j++) {
                System.out.print(" " + nf.format(gramSum[i][j]));
                if(i==0) {
                    assertEquals(2.00*Math.cos(j*Math.PI/k), gramSum[0][j], 0.001);
                }
            }
            System.out.println();
        }
        for (int j = 0; j < 2*k; j++) {
            System.out.print(" " + nf.format(hist[j].mean()));
            assertEquals(2.00*Math.cos(j*Math.PI/k), hist[j].mean(), 0.001);
        }
        System.out.println();
        for (int j = 0; j < 2*k; j++) {
            System.out.print(" " + nf.format(hist[j].stdDev()));
        }
        System.out.println();

        //quantile. use R
        final File outputZquantileFile = new File(baseDir, "gzeta_quantile" + k + ".csv");
        //final PrintWriter outputZquantile = new PrintWriter(outputZquantileFile);
        final PrintWriter outputZquantile = null;
        if(outputZquantile != null) {
            //loop over phi
            double[] f = {
                    -0.1,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1.1
            };
            outputZquantile.println(
                    "phi*6/pi,min,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,max,mean"
            );
            for (int j = 0; j < 2*k; j++) {
                outputZquantile.print(j + ", ");
                double quantile = hist[j].yForIndex(1)+hist[j].delta/2.0;
                NormalizedSpline cdfSpline = hist[j].cdfSpline();
                for(int fidx = 0; fidx < f.length; fidx++) {
                    if(fidx == 1) {
                        quantile = hist[j].findQuartile(f[fidx], epsilon);
                    } else {
                        quantile = hist[j].findQuartile(f[fidx], quantile, epsilon);
                    }
//	        		System.out.println(j + " " + quantile + " " + f[fidx]);
//	        		System.out.println(cdfSpline.eval(quantile));
                    outputZquantile.print(nf.format(quantile) + ", ");
                }
                outputZquantile.print(nf.format(hist[j].mean()) + ", ");
                outputZquantile.println();
            }
            outputZquantile.close();
        }

        //hist
        final File outputHistFile = new File(baseDir,
                "calcHist_fine" + k + ".csv");
        //		"calcHist" + k + ".csv");
        final PrintWriter outputHist = new PrintWriter(outputHistFile);
        if(outputHist != null) {
            final double norm = hist[0].sampleSize*hist[0].delta;
            if(k == 6) {
                outputHist.println(",0,30,60,90,120,150,180,210,240,270,300,330,");
            } else {
                outputHist.println(",0,15,30,45,60,75,90,105,120,135,150,165,180,195,210,225,240,255,270,285,300,315,330,345,");
            }
            for(int i = 0; i < hist[0].hist.length; i++) {
                outputHist.print(nf.format(hist[0].yForIndex(i)+hist[0].delta/2.0) + ", " );
                for (int j = 0; j < 2*k; j++) {
                    outputHist.print(nf.format((double)hist[j].hist[i]/norm) + ", ");
                }
                outputHist.println();
            }
            outputHist.close();
        }
    }

    /**
     * 398942 terms of GSeries
     * @param t0
     * @param initialPadding
     * @return
     * @throws IOException
     */
    private GSeries calculateGSeriesE12( double t0, int initialPadding) throws IOException {
        int index = (int) Math.floor(t0);
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        double begin =  t0;
        int k0 = 1, k1=398942;
        DataOutputStream out = null;
        File file = new File("out/" + Integer.toString(index) +"E12.dat");
        boolean output = false;
        if (output) {
            if (!file.getParentFile().exists()) {
                file.getParentFile().mkdirs();
            }
            try {
                OutputStream os = new FileOutputStream(file);
                BufferedOutputStream bos = new BufferedOutputStream(os);
                // create data output stream
                out = new DataOutputStream(bos);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
        long init= System.currentTimeMillis();
        double incr  = Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
        int R = 500+2*initialPadding;
        System.out.println("incr " + incr);
        begin -= initialPadding*incr;
        GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, R);
        long end = System.currentTimeMillis();
        System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
        if (output) {
            out.writeDouble(begin);
            out.writeDouble(incr);
            for (int i = 0; i < gAtBeta.gAtBeta.length; i++) {
                out.writeDouble(gAtBeta.gAtBeta[i][0]);
                out.writeDouble(gAtBeta.gAtBeta[i][1]);
            }
            out.close();
        }

        return gAtBeta;
    }

    private double[][] writeZetaPhi(int sampleIndex, PrintWriter outputZ, double[][] gramSum,
                                    Histogram[] hist) throws Exception{
        double t0 = StaticMethods.gramE12[sampleIndex][0];
        final boolean calculateCross = gramSum.length > 1;
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = StaticMethods.getSavedGSeries(t0, offset);
        //double[] oddsum = new double[k], evensum = new double[k];
        final double firstGram = Gram.gram(offset, t0 + 0.001 );
        final int N = 29999;
        long gramIndex = Gram.gramIndex(offset, firstGram);
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+firstGram)/(2*Math.PI)));
        BigDecimal tvalsi = offset.add(BigDecimal.valueOf(firstGram+ N*incr/6), Gram.mc);
        incr = Gram.gramInterval(tvalsi);
        System.out.println("****** " + (gramIndex-3945951431271L));

        double gram = firstGram-incr;
        int k = gramSum[0].length/2;
        final double[] zeta = new double[2*k];
        //final double[][] gramSum = new double[2][2*k];
        double[] saved =  new double[k];
        for (int i = 0; i < N; i++) {
            if(i>0 && i*3%N == 0){
                tvalsi = tvalsi.add(BigDecimal.valueOf(N*incr/3), Gram.mc);
                incr = Gram.gramInterval(tvalsi);
            }
            gram += incr;
            for (int j = 0; j < k; j++) {
                double t = gram + j*incr/k;
                double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset(
                        t, 4, 40, 1.6E-9, false);
                if(i%2==0){
                    // we begin with odd.
                    zeta[k+j] = gAtBeta.riemannZeta(gFromBLFI, t);
                    gramSum[0][j+k] += zeta[k+j];
                    hist[j+k].addPoint(zeta[k+j]);
                    if(calculateCross && i>0 ) {
                        gramSum[1][j+k] += saved[j]*zeta[j+k];
                    }
                } else {
                    zeta[j] = gAtBeta.riemannZeta(gFromBLFI, t);
                    gramSum[0][j] += zeta[j];
                    hist[j].addPoint(zeta[j]);

                    if(calculateCross) {
                        gramSum[1][j] += zeta[j] * zeta[j+k];
                    }
                    saved[j] = zeta[j];
                }
            }
            if(outputZ!=null && i%2==1){
                for (int j = 0; j < 2*k; j++) {
                    if(j>0){outputZ.print(", ");}
                    outputZ.print(nf.format(zeta[j]));
                }
                outputZ.println();
            }
        }
        for (int j = 0; j < 2*k; j++) {
            gramSum[0][j] *= 2.0/N;
            for (int i1 = 1; i1 < gramSum.length; i1++) {
                gramSum[i1][j]  *= 2.0/N;
            }
        }
        for (int j = 0; j < 2*k; j++) {
            assertEquals(2.00*Math.cos(j*Math.PI/k), gramSum[0][j], 0.05);
        }
        //        long actual = (gramIndex-3945951431271L)%N;
        //        assertTrue("index " + actual, actual==0 || actual==1);
//	        System.out.println(firstGram + incr*N);
        return gramSum;
    }


    @Test //@Ignore
    public void testDoubleDer() throws Exception {
        double xmid = 1;
        double delta = 0.001;
        double x = xmid;
        double f0 = x*x*x;
        x = xmid + delta;
        double f1 = x*x*x;
        x = xmid - delta;
        double f_1 = x*x*x;
        double doubleDer = (f1 + f_1 -  2*f0)/(delta*delta);
        assertEquals(6, doubleDer, 1.0E-9);
        double der = (f1 - f_1 )/(2*delta);
        assertEquals(3, der, 1.0E-6);

    }

    @Test //@Ignore
    public void testGetSavedGSeries1() throws Exception{
//		double firstZero = 243837.44036866794;
        double firstZero = 243831.456494008;
        int idx = StaticMethods.findFile(firstZero);
        double maxZeroDev = Double.MIN_VALUE;
        double maxDerDev = Double.MIN_VALUE;
        double maxMaxDev = Double.MIN_VALUE;

        GSeries gAtBeta;
        double upperforzero = 0.000001;
        double deltader = 0.000025;
        double deltamax = 5.0E-7;
        int sampleSize = 25;
        String gbetaSource = "Interpolat";

        System.out.println("** gbetaSource " + gbetaSource + " ======");
        switch (gbetaSource) {
            case "Saved":
                double t0 = StaticMethods.gramE12[idx][0];
                final BigDecimal offset = BigDecimal.valueOf(1.0E12);
                gAtBeta = StaticMethods.getSavedGSeries(t0, offset);
                break;
            case "Interpolate":
                gAtBeta = Interpolate.readGSeries();
                upperforzero = 0.5;
                deltader = 10.0;
                deltamax = 1.0;
                break;
            default:
                FixE12GSeries fixE12GSeries = new FixE12GSeries();
                fixE12GSeries.testChangeToZetaAndDer(fixE12GSeries.gAtBeta);
                gAtBeta = fixE12GSeries.gAtBeta;
                upperforzero = 0.5;
                deltader = 10.0;
                deltamax = 1.0;
                
        }

        int iMax = 0;
        int iZero = 0;
        int iDer = 0;
    
        //25 rows zero expectedDer

        double[] nextValues = CopyZeroInformation.skipUntil(Interpolate.zeroIn, firstZero);
        int i = 0;
        double z0 = 0, d0 = -1.0, extremumFromFile = -1.0;
        // cant go below 40
        final int initialPadding = 40;
    
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        
        for (i = 0; i <= sampleSize; i++) {
            double zeroPosition = nextValues[0];
            double expectedDer =  nextValues[1];
            double zeta = gAtBeta.evaluateZeta(zeroPosition, initialPadding);
            assertEquals("zero", 0.0, zeta, upperforzero);
            double der = gAtBeta.evalDer(
                    zeroPosition, initialPadding, 0.00025*gAtBeta.spacing);
            assertEquals("der",expectedDer, der, deltader);
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
                assertEquals("max", extremumFromFile, evalMax, deltamax);
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
                nextValues = CopyZeroInformation.skipUntil(Interpolate.zeroIn, nextValues[0]);
            }
        }
        System.out.println("done");
        /*
        for (double[] zeroEntry: zeroInfo) {
            System.out.println(
                Arrays.toString(zeroEntry)
            );
        }
        
         */
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
    }

    /**
     * Test calculate Z using F, G from saved file
     * @throws Exception
     */
    @Test //@Ignore
    public void testGetSavedGSeries() throws Exception{
        double t = 247.149;
        int idx = StaticMethods.findFile(t);
        final int initialPadding = 40;

        double maxDerDev = Double.MIN_VALUE;
        double maxMaxDev = Double.MIN_VALUE;
        double t0 = StaticMethods.gramE12[idx][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = StaticMethods.getSavedGSeries(t0, offset);
        {
            double zeroPosition = 243.8749480149;
            double upper = 244.0199875;
            double incr = (upper-zeroPosition)/2;
            while (zeroPosition<=upper){
                double zeta = gAtBeta.evaluateZeta(zeroPosition, 40);
                double der = gAtBeta.evaluateDer(zeroPosition, 40);
                System.out.println("t " + zeroPosition
                        + " : zeta from GSeries: " + zeta
                        + " : der from GSeries: " + der
                );
                zeroPosition += incr;

            }

        }
        //25 rows zero expectedDer
        String zerosFile = "data/gzetaE12/zerosE12.csv";
        BufferedReader zeroIn = new BufferedReader(new FileReader(zerosFile));
        String in = zeroIn.readLine();
        for (int i = 0; i < 1; i++) {
            in = zeroIn.readLine();
        }

        int i = 0;
        double z0 = 0, d0 = -1.0, extremumFromFile = -1.0;
        while (in != null) {
            String[] dataFromFile = in.split(",\\s+");
            double zeroPosition = Double.parseDouble(dataFromFile[0]);
            double expectedDer = Double.parseDouble(dataFromFile[1]);
            double zeta = gAtBeta.evaluateZeta(zeroPosition, 40);
            assertEquals(0.0, zeta, 0.000001);
            double der = gAtBeta.evaluateDer(zeroPosition, 40);
            assertEquals(expectedDer, der, 0.00006);
            System.out.println("** i " + ++i);
            System.out.println("{" + in + "},");
            System.out.println("zeroPosition " + zeroPosition + " : eval from GSeries: " + zeta);
            System.out.println(
                    "expectedDer  "
                            + expectedDer
                            + " : eval from GSeries: " + der
                            + " diff " + Math.abs(expectedDer-der)
            );
            if (i>1) {
                Poly4 poly = new Poly4(z0, zeroPosition, d0, expectedDer,
                        extremumFromFile);
                double positionMax = poly.getPositionMax();
                double evalMax = gAtBeta.evaluateZeta(positionMax, 40);
                double maxDev = extremumFromFile - evalMax;
                double[] maxder = gAtBeta.doubleDer(positionMax, initialPadding,
                        evalMax, 0.0005*gAtBeta.spacing);
                System.out.println(
                        "positionMax " + positionMax
                                + ", eval " + evalMax
                                + " read " + extremumFromFile
                                + " diff(Max) " + maxDev
                );
                for (int j = 0; j < 2; j++) {
                    if(Math.abs(maxder[0]) > 0.0001 &&  Math.abs(maxDev) > 0.00000001) {
                        positionMax -= maxder[0] /maxder[1];
                        evalMax = gAtBeta.evaluateZeta(positionMax, initialPadding);
                        maxDev = extremumFromFile - evalMax;
                        System.out.println(
                                "positionMax " + positionMax
                                        + ", eval " + evalMax
                                        + " read " + extremumFromFile
                                        + " diff(Max) " + maxDev
                        );
                        maxder = gAtBeta.doubleDer(positionMax, initialPadding,
                                evalMax, 0.0005*gAtBeta.spacing);
                        System.out.println(
                                "positionMax der " + Arrays.toString(maxder)
                        );
                    }
                }

                if(Math.abs(maxDev) > maxMaxDev){
                    maxMaxDev = Math.abs(maxDev);
                }
                assertEquals(extremumFromFile, evalMax, 5.0E-7);
            }
            z0 = zeroPosition;
            d0 = expectedDer;
            extremumFromFile = Double.parseDouble(dataFromFile[2]);
            in = zeroIn.readLine();
        }
        zeroIn.close();
        System.out.println("done");
        System.out.println(
                "maxMaxDev  " + maxMaxDev
        );
    }

    @Test //@Ignore
    public void test1E12Zeros() throws Exception {
		/*
		compare calculated with saved
		 */
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        double t0 = 244.021159171564;
        System.out.println(nf.format(t0) + ", " + Gram.gramIndex(offset, t0));
        double zeroFromZerosFile = 243.8749480149;
        double expectedDerFromZerosFile = 17.619276585379914;
        double oldZeroFromZerosFile = 243.8749480149;
        double oldExpectedDerFromZerosFile = 17.619276585379914;
        double extremumFromFile = 1.9266754104451154;
        int currentIndex = 40;
        final int initialPadding = currentIndex;
        // calculate GSeries beginning at t0
        GSeries gAtBeta = calculateGSeriesE12(t0, initialPadding);

        double[] gAtGram = gAtBeta.gAtBeta[currentIndex];
        Pair<double[], Double> fAndZGram = gAtBeta.fAndZ(gAtGram, t0);
        System.out.println("At " + t0 + ", " + fAndZGram.getValue() + ", " + gAtBeta.spacing);
        assertEquals(1.92649807303971, fAndZGram.getValue(), 0.0000015);

        double correction = gAtBeta.correctionAtT(zeroFromZerosFile);
        String zetaFile = "data/zetaE12.csv";
        BufferedReader zetaIn =
                new BufferedReader(new FileReader(zetaFile));
        String Zinput = zetaIn.readLine();
        for (int i = 0; i < 2; i++) {
            Zinput = zetaIn.readLine();
        }
        String zerosFile = "data/gzetaE12/zerosE12.csv";
        BufferedReader zeroIn = new BufferedReader(new FileReader(zerosFile));
        String inFromZerosFile = zeroIn.readLine();
        // this line should be empty.
        for (int i = 0; i < 1; i++) {
            inFromZerosFile = zeroIn.readLine();
            //first saved zero
        }

        int gramIndex = Integer.MAX_VALUE;
        while (inFromZerosFile != null) {
            String[] splitZerosLine = inFromZerosFile.split(",\\s+");
            zeroFromZerosFile = Double.parseDouble(splitZerosLine[0]);
            expectedDerFromZerosFile = Double.parseDouble(splitZerosLine[1]);
            {
                Poly4 poly = new Poly4(oldZeroFromZerosFile, zeroFromZerosFile,
                        oldExpectedDerFromZerosFile, expectedDerFromZerosFile,
                        extremumFromFile);
                double positionMax = poly.getPositionMax();
                double evalMax = gAtBeta.evaluateZeta(positionMax, 40);
                System.out.println(
                        "positionMax " + positionMax
                                + ", eval " + evalMax
                                + " read " + extremumFromFile
                                + " diff(Max) " + (extremumFromFile-evalMax)
                );
                assertEquals(extremumFromFile, evalMax, 0.013);
            }
            oldZeroFromZerosFile = zeroFromZerosFile;
            oldExpectedDerFromZerosFile = expectedDerFromZerosFile;
            extremumFromFile = Double.parseDouble(splitZerosLine[2]);

            while(zeroFromZerosFile>t0) {
                String[] parsed = Zinput.split(",");
                double zetaSaved = Double.parseDouble(parsed[1]);
                gramIndex = Integer.parseInt(parsed[0]);
                assertEquals(zetaSaved, fAndZGram.getValue(), 0.000001);
                Zinput = zetaIn.readLine();

                System.out.println("**Gram " + gramIndex + " " + t0 + " " + fAndZGram.getValue());
                System.out.println("*F at Gram " + Arrays.toString(fAndZGram.getKey()));
                t0 += 2 * gAtBeta.spacing;
                currentIndex += 2;
                gAtGram = gAtBeta.gAtBeta[currentIndex];
                fAndZGram = gAtBeta.fAndZ(gAtGram, t0);
            }

            double[] gFromBLFI0 = gAtBeta.diagnosticBLFISumWithOffset(zeroFromZerosFile, 4,
                    initialPadding, 1.6E-9, false);
            double theta = gAtBeta.theta(zeroFromZerosFile);
            double cos = Math.cos(theta);
            double sin = Math.sin(theta);

            Pair<double[], Double> fAndZero = gAtBeta.fAndZ(gFromBLFI0, zeroFromZerosFile);
            double zeta = fAndZero.getValue();
            double[] f = fAndZero.getKey();
            //
            System.out.println("zeroFromZerosFile " + zeroFromZerosFile + " zeta " + nf.format(zeta ) + " cf 0.0");
            System.out.println("*F at zeroFromZerosFile" +   Arrays.toString(fAndZero.getKey()));
            assertEquals(0.0, zeta, 0.000001);

            Pair<double[], Double> fZprime = gAtBeta.der(zeroFromZerosFile, initialPadding);
            double der = fZprime.getValue();
            double[] fprime = fZprime.getKey();
            double zPrime = 2*(cos*fprime[0] + sin*fprime[1])
                    + 2*gAtBeta.lnsqrtArg1*(-sin*f[0] + cos*f[1])
                    ;
            assertEquals(expectedDerFromZerosFile, zPrime, 0.00003);
            assertEquals(expectedDerFromZerosFile, der, 0.00003);
            // read zero from line
            inFromZerosFile = zeroIn.readLine();
//			if (gramIndex>1) {
//				break;
//			}
        }

        zeroIn.close();
        zetaIn.close();
    }

    /**
     * Test method for fseries at gram points 1.0E12.
     * Exact calculation, compare with riemann.InterpolateTest.testReadItems()
     * @throws FileNotFoundException
     */
    @Test @Ignore
    public void test1E12F() throws FileNotFoundException {
        int k0 = 1, k1=398942;
        Gram.initLogVals(k1);
        int R = 20;
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        //gram [GramInfo [grampt=1000000000243.77756012466052947405878015472510, idx=3945951431271], 
        //      GramInfo [grampt=1000000000244.02115917156451839965694310614387, idx=3945951431272]]
        //0.24359904690398881 tincr 0.24359904690399015
        double incr  = 0.24359904590398668;
        double begin = 243.77756012466052947405878015472510 + 2*incr;
        BigDecimal tval = new BigDecimal(begin, Gram.mc).add(
                offset, Gram.mc);
        double[][] fAtBeta = GSeries.fSeries(k0, k1, incr/2, 2*R, tval);

        BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
        double predictedSqrtArg1 = sqrtArg1.doubleValue();
        double cos = -1;
        double correction = GSeries.correction( predictedSqrtArg1);

        File file = new File("out/gSeriesE12/fseries.csv");
        PrintStream out = new PrintStream(file);
        //PrintStream out = System.out;
        out.println("n-3945951431270L,  f.real, f.im, zeta");
        for (int i = 0; i < R; i++) {
            double rotatedSum = 2*( cos*fAtBeta[2*i][0]);
            double zeta = rotatedSum + correction;
            cos = -cos;
            out.println((i+3) + ", " + nf.format(fAtBeta[2*i][0]) + ", "
                    + nf.format(fAtBeta[2*i][1]) + ", " + nf.format(zeta));
            out.println((i+3) + ", " + nf.format(fAtBeta[2*i+1][0]) +
                    ", " + nf.format(fAtBeta[2*i+1][1]));
        }
        out.close();
    }

    private void storeGE12(double zero, double t0) throws IOException, FileNotFoundException {
        int index = (int) Math.floor(t0);
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        double begin = Gram.gram(offset, t0 );
        int k0 = 1, k1=398942;
        DataOutputStream out = null;
        File file = new File("out/gSeriesE12/" + Integer.toString(index) +".dat");
        boolean output = true;
        if (output) {
            if (!file.getParentFile().exists()) {
                file.getParentFile().mkdirs();
            }
            try {
                OutputStream os = new FileOutputStream(file);
                BufferedOutputStream bos = new BufferedOutputStream(os);
                // create data output stream
                out = new DataOutputStream(bos);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
        long init= System.currentTimeMillis();
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
        final int initialPadding = 40;
        int R = 30000+2*initialPadding;
        System.out.println(incr);
        begin -= initialPadding*incr;
        GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, R);
        long end = System.currentTimeMillis();
        System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
        if (output) {
            out.writeDouble(begin);
            out.writeDouble(incr);
            for (int i = 0; i < gAtBeta.gAtBeta.length; i++) {
                out.writeDouble(gAtBeta.gAtBeta[i][0]);
                out.writeDouble(gAtBeta.gAtBeta[i][1]);
            }
            out.close();
        }
        System.out.println(gAtBeta.riemannZeta(gAtBeta.gAtBeta[initialPadding], begin));
        //g  : [-0.33143958775035764, 0.0733285174786178] 1287.5146091794
        double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, initialPadding, 1.6E-9, false);
        double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
        System.out.println("** g  : " + Arrays.toString(gFromBLFI) + " zeta " + zeta);
        assertTrue(Math.abs(zeta) < 0.000001);
    }

    private int testCorrelationE12(int sampleIndex, PrintWriter out, double[][] prod) throws Exception{
        double t0 = StaticMethods.gramE12[sampleIndex][0];
        final BigDecimal offset = BigDecimal.valueOf(1.0E12);
        GSeries gAtBeta = StaticMethods.getSavedGSeries(t0, offset);
        int k = prod.length/2;
        double[] oddsum = new double[k], evensum = new double[k];
        final double[] zeta = new double[2*k];
        final double firstGram = Gram.gram(offset, t0 + 0.001 );
        final int N = 30000;
        long gramIndex = Gram.gramIndex(offset, firstGram);
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+firstGram)/(2*Math.PI)));
        BigDecimal tvalsi = offset.add(BigDecimal.valueOf(firstGram+ N*incr/6), Gram.mc);
        incr = Gram.gramInterval(tvalsi);
        System.out.println("****** " + (gramIndex-3945951431271L));

        double gram = firstGram;
        int size = 0;
        for (int i = 0; i < N-1; i++) {
            if(i>0 && i*3%N == 0){
                tvalsi = tvalsi.add(BigDecimal.valueOf(N*incr/3), Gram.mc);
                incr = Gram.gramInterval(tvalsi);
            }
            gram += incr;
            for (int j = 0; j < k; j++) {
                double t = gram + j*incr/k;
                double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset(
                        t, 4, 40, 1.6E-9, false);
                if(i%2==0){
                    zeta[j] = gAtBeta.riemannZeta(gFromBLFI, t);
                    evensum[j] += zeta[j];
                } else {
                    zeta[k+j] = gAtBeta.riemannZeta(gFromBLFI, t);
                    oddsum[j] += zeta[k+j];
                }
            }
            if(i%2==1){
                for (int i1 = 0; i1 < 2*k; i1++) {
                    for (int j = i1; j < 2*k; j++) {
                        prod[i1][j] += zeta[i1]*zeta[j];
                    }
                }
                size++;
            }
        }
        for (int j = 0; j < k; j++) {
            assertEquals(-2.00*Math.cos(j*Math.PI/k), oddsum[j]/size, 0.05);
            assertEquals(2.00*Math.cos(j*Math.PI/k), evensum[j]/size, 0.05);
        }
        assertEquals(0, (gramIndex-3945951431271L)%N);
        return size;
    }

    /**
     * Compare zeta from F-series with zeta from Riemann evaluation.
     * Z(t) = Real(exp(−i*theta(t))F(1,floor(tau); t)) + R(t)
     * dtheta/dt = ln(t/(2pi))/2;
     */
    @Test @Ignore
    public void testZeroLargeOffset() {
        double[][] fAtBeta = null;
        double[] begin = {243.8749480149, 1436233.106281030331450810};
        int k0 = 1, k1=0;
        int R = 2;
        double lnsqrtArg1 = 0;
        double basetheta = 0;
        double dsqrtArg1 = 0;
        double tbase = 0;
        double basesqrtArg1 = 0;
        long offset = (long) 1.0E12;
        for (int i = 0; i < begin.length; i++) {
            double tincr =  (begin[i]-begin[0]) ;
            BigDecimal tval = new BigDecimal(begin[i], Gram.mc).add(
                    BigDecimal.valueOf(offset), Gram.mc);
            double predictedSqrtArg1 = 0;
            double theta = 0;
            if(i == 0 ){
                dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
                tbase = tval.doubleValue();
                BigDecimal t2 = tval.divide(Gram.bdTWO);
                BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
                basesqrtArg1 = sqrtArg1.doubleValue();
                k1 = (int)basesqrtArg1;
                BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, Gram.mc);
                lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
                fAtBeta = GSeries.fSeries(k0, k1, begin[1]-begin[0], R, tval);
                //theta should be 2.819633653651107
                theta = tval.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
                        .subtract(Gram.pi8, Gram.mc).remainder(Gram.pi_2).doubleValue();
                basetheta = theta;
                predictedSqrtArg1 = basesqrtArg1 ;
                //0.24359904690398668
                GramInfo[] gram = GramInfo.RZGram(3945951431271L, 2, Gram.mc);
                System.out.println("gram " + Arrays.toString(gram));
            } else {
                theta = (basetheta + lnsqrtArg1*tincr
                        +tincr*tincr/(4*tbase))%(2*Math.PI);
                predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;

                BigDecimal tBase = new BigDecimal(begin[i], Gram.mc).add(
                        BigDecimal.valueOf(offset), Gram.mc);
                BigDecimal alphaBD = (Gram.log(k0).add(Gram.log(k1))).divide(Gram.bdTWO, Gram.mc);
                double argalphaBase = tBase.multiply(alphaBD, Gram.mc).remainder(Gram.pi_2).doubleValue();
                double[] g = new double[2];
                g[0] = Math.cos(argalphaBase)*fAtBeta[i][0] + Math.sin(argalphaBase)*fAtBeta[i][1];
                g[1] = Math.cos(argalphaBase)*fAtBeta[i][1] - Math.sin(argalphaBase)*fAtBeta[i][0];
                System.out.println("g  : " + Arrays.toString(g) + " argalphaBase " + argalphaBase);

            }
            double rotatedSum = 2*( Math.cos(theta)*fAtBeta[i][0]+Math.sin(theta)*fAtBeta[i][1]);
            double correction = GSeries.correction( predictedSqrtArg1);
            double zeta = rotatedSum + correction;
            System.out.println("f  : " + Arrays.toString(fAtBeta[i])
                    + " theta " + theta + " rotatedSum " + rotatedSum
                    + " zeta " + zeta);
            double zetaFromRiemann = Riemann.riemann(begin[i], offset);
            System.out.println("zetaFromRiemann " + zetaFromRiemann);
            System.out.println("sqrtArg1[i].doubleValue() " + predictedSqrtArg1 + " correction " + correction );
            assertTrue(i + " ", Math.abs(zeta)  < 5.0E-7);
            assertTrue(i + " ", Math.abs(zeta-zetaFromRiemann)  < 5.0E-7);
        }
    }


    /**
     * Test method for small values of t (which don't need BigDecimal).
     * Change G Value At Index
     */
    @Test
    public void testGSeriesSmallTIncrementGSeries() {
        int k0 = 10, k1=100;
        int N = 30;
        long init= System.currentTimeMillis();
        GSeries gSeries = new GSeries(k0, k1, 5, 5+N-1);
        long end = System.currentTimeMillis();
        System.out.println("calc for " + N + ": " + (end - init) + "ms");
        int minIndex = 5;
        double t0 = (minIndex+N/2+0.5)*gSeries.spacing;
        double[] gFromBLFI = gSeries.blfiSumWithOffsetSmallT( t0, 2);
        int midIdx = gSeries.midIdx;
        System.out.println("midIdx " + midIdx);
        double[] directEval = gSeries.gSeriesForSmallT(t0);
        assertTrue(Math.abs(gFromBLFI[0] - directEval[0])
                + Math.abs(gFromBLFI[1] - directEval[1]) < 0.005);
        System.out.println(" t0 " +  t0 + " gFromBLFI 0 " + gFromBLFI[0] + ", 1 " + gFromBLFI[1]
                + ": directEval" + Arrays.toString(directEval));
        gSeries.incrementGValueAtIndex(midIdx, new double[]{100, 1000});
        gFromBLFI = gSeries.blfiSumWithOffsetSmallT( t0, 2);
        double factorAtIndex = gSeries.factorAtIndex(midIdx, t0, 2);
        System.out.println(" t0 " +  t0 + " gFromBLFI 0 " + gFromBLFI[0] + ", 1 " + gFromBLFI[1]
                + ", factorAtIndex "
                + factorAtIndex);
        assertTrue(Math.abs(gFromBLFI[0] - directEval[0]-100*factorAtIndex)
                + Math.abs(gFromBLFI[1] - directEval[1]-1000*factorAtIndex) < 0.005);
    }

}
