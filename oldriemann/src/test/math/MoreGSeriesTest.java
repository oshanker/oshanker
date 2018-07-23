package math;

import static org.junit.Assert.*;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Gram;
import riemann.Rosser;

public class MoreGSeriesTest {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(8);
        nf.setMaximumFractionDigits(8);
        nf.setGroupingUsed(false);
    }
    
    @Test @Ignore
    public void testSymmetryRelations() throws Exception{
        //check the symmetry and antisymmetry relations from the output of distributions
        File file = new File("out/kurtosis6.txt");
        BufferedReader zeroIn = new BufferedReader(new FileReader(file));
        int k = 12;
        double[][] vals = new double[k][3];
        for (int i = 0; i < k; i++) {
            String in = zeroIn.readLine();
            String[] line = in.split("\\s+");
            for (int j = 0; j < vals[i].length; j++) {
                vals[i][j]= Double.parseDouble(line[j+1].trim());
            }
        }
        double skewSum = 0;
        for (int i = 0; i < k; i++) {
            double skew = vals[i][1];
            skewSum += skew;
        }
        skewSum /= k;
        System.out.println("skewSum " + skewSum);
        double skewCoefficient = (vals[0][1]-skewSum);
        for (int i = 0; i < k; i++) {
            double skew = vals[i][1];
            System.out.println(
                    (skew-skewSum)+ ", " + skewCoefficient*Math.cos(i*2.0*Math.PI/k)
                    );
        }
       
        zeroIn.close();
    }

    @Test
    public void testInterpolate() throws Exception{
        int N = 199;
        GSeries gSeries = createGseries(N);
        final int initialPadding = 40;
        //[261.0016582858428, 28.452144305679546, 2.3693797179877887], 261.07309238484356, 
        //[261.31522681873514, -6.883771986248166, 0.08672183144103376]
        //        1.7592847984372848 ** 0.03495629927494104
//        double zeta = evaluateZeta(261.07309238484356, initialPadding, gSeries);
//        System.out.println( " zeta at Gram " + zeta + " cf 1.7592847984372848" );
//        double[] zero = {261.0016582858428, 261.31522681873514};
//        double[] expectedDer = {28.452144305679546, -6.883771986248166};
        
        /*
[109.9434127500521, 207.28544365034014, 8.283292860041835], 109.99801991618585, [110.10427375389713, -61.091725512779625, 0.8755383742860865]
7.852770303334955 ** 

[115.21645409737458, -7.282653909337562, 0.8259045475747673], 115.31471904908707, [115.35911882837084, 17.960412142999786, 0.38874142446558396]
-0.7183317846044628 ** 
         */
        double zeta = evaluateZeta(109.99801991618585, initialPadding, gSeries);
        System.out.println( " zeta at Gram " + zeta + " cf 7.852770303334955" );
        double[] zero = {109.9434127500521, 110.10427375389713, 115.21645409737458, 115.35911882837084};
        double[] expectedDer = {207.28544365034014, -61.091725512779625, -7.282653909337562, 17.960412142999786};
        for (int i = 0; i < zero.length; i++) {
            validateZero(zero[i], expectedDer[i], initialPadding, gSeries,false);
        }

//        File file = new File("out/gzetaE28/gzeta6.csv");
//        if (!file.getParentFile().exists()) {
//            file.getParentFile().mkdirs();
//        }
//        PrintWriter out = new PrintWriter(file);
//        writeZetaPhi(N, out, initialPadding);
//        out.close();
    }

    public static GSeries createGseries(int N) throws IOException, FileNotFoundException {
        Rosser.readConfig("data/RosserConfig.txt");
        int correction = 0;
        if(Rosser.configParams.containsKey("correction")){
            correction = Rosser.getParamInt("correction");
        }
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        int noffset = Rosser.getParamInt("noffset");
        final BigDecimal offset = new BigDecimal(Rosser.getParam("bdoffset"));
        double begin= baseLimit + (noffset-correction)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        double zetaCorrection = GSeries.correction( gSeries.basesqrtArg1);
        System.out.println( gSeries.begin + ", " + zetaCorrection);
        File reFile = new File(Rosser.getParam("conjecturesOutFile").replace(
                "stats", "reF_Gram_"));
        File imFile = new File(Rosser.getParam("conjecturesOutFile").replace(
                "stats", "imF_Gram_"));
        BufferedReader reFReader = new BufferedReader(new FileReader(reFile));
        BufferedReader imFReader = new BufferedReader(new FileReader(imFile));
        reFReader.readLine();
        imFReader.readLine();
        double[][] fAtBeta = new double[N-1][2];
        for (int i = 0; i < fAtBeta.length; i++) {
            String[] reF = reFReader.readLine().split(",");
            fAtBeta[i][0] = Double.parseDouble(reF[1].trim());
            String[] imF = imFReader.readLine().split(",");
            fAtBeta[i][1] = Double.parseDouble(imF[1].trim());
        }
        gSeries.rotateFtoG(fAtBeta);
        gSeries.setgAtBeta(fAtBeta);
        imFReader.close();
        reFReader.close();
        return gSeries;
    }

    private void writeZetaPhi(int N, PrintWriter out, int initialPadding) throws Exception{
        GSeries gAtBeta = createGseries(N);
        double[] oddsum = {0, 0, 0, 0, 0, 0}, evensum = {0, 0, 0, 0, 0, 0};
        int k = oddsum.length;
        final double[] zeta = new double[2*k];
        double incr  = gAtBeta.spacing;
        final double firstGram = gAtBeta.begin + initialPadding*incr;
        N -= 2*initialPadding;
        double gram = firstGram-incr;
        for (int i = 0; i < N; i++) {
            gram += incr;
            for (int j = 0; j < k; j++) {
                double t = gram + j*incr/k;
                double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( 
                        t, 4, initialPadding, 1.6E-9, false);
                if(i%2==1){
                    zeta[k+j] = gAtBeta.riemannZeta(gFromBLFI, t);
                    oddsum[j] += zeta[k+j];
                } else {
                    zeta[j] = gAtBeta.riemannZeta(gFromBLFI, t);
                    evensum[j] += zeta[j];
                }
            }
            if(i%2==1){
                for (int j = 0; j < 2*k; j++) {
                    if(j>0){out.print(", ");}
                    out.print(nf.format(zeta[j]));
                }
                out.println();
            }
        }
        for (int i = 0; i < oddsum.length; i++) {
            oddsum[i] *= 2.0/N;
            evensum[i] *= 2.0/N;
        }
        System.out.println(Arrays.toString(oddsum));
        System.out.println(Arrays.toString(evensum));
        for (int j = 0; j < k; j++) {
            assertEquals(-2.00*Math.cos(j*Math.PI/k), oddsum[j], 0.1);
            assertEquals(2.00*Math.cos(j*Math.PI/k), evensum[j], 0.1);
        }
    }
    
    @Test @Ignore 
    public void test1E12() throws Exception{
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        double t0 = 244.021159171564;
        long gramIndex = Gram.gramIndex(offset, t0);
        System.out.println(nf.format(t0) + ", " + gramIndex);
        double zero = 244.920599505825;
        double expectedDer = 23.85164367971759;
        final int initialPadding = 40;
        GSeries gAtBeta = testE12(t0, initialPadding);
        double riemannZeta = gAtBeta.riemannZeta( gAtBeta.gAtBeta[40], t0);
        System.out.println( t0 + ", " +riemannZeta);
        assertEquals(1.92649807303971, riemannZeta, 0.0000015);
        validateZero(zero, expectedDer, initialPadding, gAtBeta,true);
    }

    private GSeries testE12( double t0, int initialPadding) throws IOException, FileNotFoundException {
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

    public static void validateZero(double zero, double expectedDer, final int initialPadding, 
            GSeries gAtBeta, boolean checkAssert) {
        double zeta = evaluateZeta(zero, initialPadding, gAtBeta);
        System.out.println( " zeta " + zeta + " cf 0.0" );
        if(checkAssert){
        assertTrue(Math.abs(zeta) < 0.000001);
        }
        
        double delta = 0.01*gAtBeta.spacing;
        double zetaplus = evaluateZeta(zero+delta, initialPadding, gAtBeta);
        double zetaminus = evaluateZeta(zero-delta, initialPadding, gAtBeta);
        //244.92059950582586, 23.85164367971759, 1.0396728565623998
        double der = (zetaplus-zetaminus)/(2*delta);
        System.out.println("der " + der + " cf " + expectedDer);
        if(checkAssert){
        assertEquals(expectedDer,der, 0.001);
        }
    }

    public static double evaluateZeta(double zero, final int initialPadding, GSeries gAtBeta) {
        double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, 
                initialPadding, 1.6E-9, false);
        System.out.println("params " + zero + ", " + gAtBeta.midIdx 
                + Arrays.toString(gAtBeta.gAtBeta[gAtBeta.midIdx]));
        double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
        return zeta;
    }
}
