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
import riemann.Interpolate;

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
    public void testFix() throws Exception{
        GSeries gSeries = Interpolate.readGSeries();
        System.out.println("gSeries.begin " + gSeries.begin);
        final int initialPadding = 40;
        double zero = 100798.08697164342;
        double expectedDer = 83.55187028339371;
        double[] tmin = {
                100798.11938493361, 1.2667769416536376,
                };
        NumberFormat nf1 = NumberFormat.getInstance();

        nf1.setMinimumFractionDigits(2);
        nf1.setMaximumFractionDigits(2);
        nf1.setGroupingUsed(false);

        for (int i = 0; i < tmin.length; i++) {
            double d = tmin[i++];
            double expected = tmin[i];
            double min = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            System.out.println(nf.format(d) + ", " + nf.format(min) + ", " 
            + nf.format(expected)
            + ", " + nf1.format((min-expected)*100.0/expected)  
            + ", " + nf1.format((min-expected)));
        }
        
        int midIdx1 = -1;
        double[] g0incr = FixGSeries.evalGSeriesIncrement(gSeries, midIdx1, initialPadding, 
                new double[]{zero, expectedDer, tmin[0], tmin[1]});
        midIdx1 = gSeries.midIdx;
        gSeries.incrementGValueAtIndex(midIdx1, g0incr);
        System.out.println(Interpolate.evaluateZeta(zero, initialPadding-1, gSeries)
               + ", "  + Interpolate.evaluateDer(zero, initialPadding-1, gSeries));

        for (int i = 0; i < tmin.length; i++) {
            double d = tmin[i++];
            double expected = tmin[i];
            double min = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            System.out.println(nf.format(d) + ", " + nf.format(min) + ", " 
            + nf.format(expected)
            + ", " + nf1.format((min-expected)*100.0/expected)  
            + ", " + nf1.format((min-expected)));
        }
    }

    @Test @Ignore
    public void testFindMin() throws Exception{
        GSeries gSeries = Interpolate.readGSeries();
        final int initialPadding = 40;
        double tmin = 100798.11938493361;
        double expectedMin = 1.2667769416536376;
        double xmin = FixGSeries.xmin(gSeries, initialPadding, tmin);
        System.out.println("xmin " + xmin);
        double evalAtMin = Interpolate.evaluateZeta(xmin, initialPadding, gSeries);
        System.out.println(evalAtMin + ", cf " + expectedMin);
    }

    @Test @Ignore
    public void testMin() throws Exception{
        GSeries gSeries = Interpolate.readGSeries();
        System.out.println("gSeries.begin " + gSeries.begin);
        final int initialPadding = 40;
        double[] tmin = {
                480.8478824361364, -0.13155469450002938, 
                481.0382017785302, -1.7233926754713067, 
                100415.55752941193, -1.543823256271554, 
                100797.9832060741, -13.139079349881037,
                100798.11938493361, 1.2667769416536376,
                };
        NumberFormat nf1 = NumberFormat.getInstance();

        nf1.setMinimumFractionDigits(2);
        nf1.setMaximumFractionDigits(2);
        nf1.setGroupingUsed(false);

        for (int i = 0; i < tmin.length; i++) {
            double d = tmin[i++];
            double expected = tmin[i];
            double min = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            System.out.println(nf.format(d) + ", " + nf.format(min) + ", " 
            + nf.format(expected)
            + ", " + nf1.format((min-expected)*100.0/expected)  
            + ", " + nf1.format((min-expected)));
        }
    }
    
    @Test @Ignore
    public void testInterpolate() throws Exception{
        GSeries gSeries = Interpolate.readGSeries();
        final int initialPadding = 40;
        
        double[] zero = {100415.50500735927, 100415.61036506912, 
                100797.8878505715, 100798.08697164342,  };
        double[] expectedDer = {-46.06567120662985, 45.21334158268663, 
                -152.8048262150694, 83.55187028339371, };
        for (int i = 0; i < zero.length; i++) {
            Interpolate.validateZero(zero[i], expectedDer[i], initialPadding, gSeries,false);
        }

//        File file = new File("out/gzetaE28/gzeta6.csv");
//        if (!file.getParentFile().exists()) {
//            file.getParentFile().mkdirs();
//        }
//        PrintWriter out = new PrintWriter(file);
//        writeZetaPhi( out, initialPadding);
//        out.close();
    }

    private void writeZetaPhi( PrintWriter out, int initialPadding) throws Exception{
        GSeries gAtBeta = Interpolate.readGSeries();
        double[] oddsum = {0, 0, 0, 0, 0, 0}, evensum = {0, 0, 0, 0, 0, 0};
        int k = oddsum.length;
        final double[] zeta = new double[2*k];
        double incr  = gAtBeta.spacing;
        final double firstGram = gAtBeta.begin + initialPadding*incr;
        int N = gAtBeta.gAtBeta.length - 2*initialPadding;
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
        Interpolate.validateZero(zero, expectedDer, initialPadding, gAtBeta,true);
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
}
