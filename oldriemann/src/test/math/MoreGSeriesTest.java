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
    public void testMin() throws Exception{
        GSeries gSeries = Interpolate.readGSeries();
        final int initialPadding = 40;
        double[] tmin = {
                476.827241802265, 5.408071801391483, 
                476.9268713015278, -0.06471682573552684, 
                476.994363522793, 1.0351545777476543, 
                477.0821162228008, -0.8014217096823656, 
                477.18226212801414, 2.042944683870103, 
                477.2855533567503, -1.1214228492178184, 
                477.3888262295176, 1.7887817315646544, 
                477.4933548969591, -0.676254120819402, 
                477.5700053309393, 0.2411835034741608, 
                477.65203328367323, -0.8922806614823013, 
                477.76496430224967, 2.1582291336904778, 
                477.8666581230465, -0.5291699048804024, 
                477.96747134986134, 1.9319290719186766, 
                478.0601196766958, -0.1144748520238274, 
                478.1403952953824, 1.5118597103733398, 
                478.2629298628662, -4.15528629473845, 
                478.37605052724103, 1.0915774559440978, 
                478.4733728729896, -2.0207142268183143, 
                478.5639873577581, 0.5035945924176821, 
                478.67379272735695, -3.6679096422807893, 
                478.7845949945954, 0.46141426854030065, 
                478.87166174584445, -1.4918414129224056, 
                478.9701961284937, 1.092208541701172, 
                479.0842315831668, -2.7022899740961157, 
                479.20381930012854, 0.7696647070983229, 
                479.2817208357934, -0.05847695495237579, 
                479.36328406455084, 0.6904227596876762, 
                479.4555501553353, -0.3363574372285997, 
                479.56119058819473, 1.3303746648686758, 
                479.6601525656017, -0.14119436474388558, 
                479.7346538394949, 0.622890162575809, 
                479.83918472522515, -1.7152694353098337, 
                479.9688822642727, 3.2578539090179888, 
                480.08651493961855, -0.3646876657261998, 
                480.1470954242631, 0.13891792296173466, 
                480.25414036179365, -2.177834683804443, 
                480.37097452953486, 0.6894826155240571, 
                480.46987232063907, -0.8704364507106825, 
                480.54961796935777, 0.014080938747092486, 
                480.6420319822211, -1.3743524906950464, 
                480.76121009188597, 0.9712054749624994, 
                480.8478824361364, -0.13155469450002938, 
                480.92785544500765, 0.7883292284262932, 
                481.0382017785302, -1.7233926754713067, 
                };
        /*
[109.9434127500521, 207.28544365034014, 8.283292860041835], 
110.09833499416513, [110.10427375389713, -61.091725512779625, 0.8755383742860865]
gram 0.21737417505040368, 110.09833499416513 (99)
positionMax 110.13482549953413, -0.8755383742860865
mid -0.3721158150659786, 110.14849253315477 (99)
positionMax 110.21670545965252, 0.5730931457175008
         */
        NumberFormat nf1 = NumberFormat.getInstance();

            nf1.setMinimumFractionDigits(2);
            nf1.setMaximumFractionDigits(2);
            nf1.setGroupingUsed(false);
        
        for (int i = 0; i < tmin.length; i++) {
            double d = tmin[i++];
            double expected = tmin[i];
            double min = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            System.out.println(nf.format(d) + ", " + nf.format(min) + ", " + nf.format(expected)
            + ", " + nf1.format((min-expected)*100.0/expected)  + ", " + nf1.format((min-expected)));
        }
    }
    
    @Test @Ignore
    public void testInterpolate() throws Exception{
        //GSeries gSeries = Interpolate.createGseries(N);
        GSeries gSeries = Interpolate.readGSeries();
        final int initialPadding = 40;
        //[261.0016582858428, 28.452144305679546, 2.3693797179877887], 261.07309238484356, 
        //[261.31522681873514, -6.883771986248166, 0.08672183144103376]
        //        1.7592847984372848 ** 0.03495629927494104
        
        /*
[109.9434127500521, 207.28544365034014, 8.283292860041835], 109.99801991618585, [110.10427375389713, -61.091725512779625, 0.8755383742860865]
7.852770303334955 ** 

[115.21645409737458, -7.282653909337562, 0.8259045475747673], 115.31471904908707, [115.35911882837084, 17.960412142999786, 0.38874142446558396]
-0.7183317846044628 ** 
         */
        double zeta = Interpolate.evaluateZeta(109.99801991618585, initialPadding, gSeries);
        System.out.println( " zeta at Gram " + zeta + " cf 7.852770303334955" );
        double[] zero = {109.9434127500521, 110.10427375389713, 115.21645409737458, 115.35911882837084};
        double[] expectedDer = {207.28544365034014, -61.091725512779625, -7.282653909337562, 17.960412142999786};
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
