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

    @Test
    public void testInterpolate() throws Exception{
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
                "stats", "reFGram"));
        File imFile = new File(Rosser.getParam("conjecturesOutFile").replace(
                "stats", "imFGram"));
        BufferedReader reFReader = new BufferedReader(new FileReader(reFile));
        BufferedReader imFReader = new BufferedReader(new FileReader(imFile));
        reFReader.readLine();
        imFReader.readLine();
        int N = 2000;
        double[][] fAtBeta = new double[N-1][2];
        for (int i = 0; i < fAtBeta.length; i++) {
            String[] reF = reFReader.readLine().split(",");
            fAtBeta[i][0] = Double.parseDouble(reF[1].trim());
            String[] imF = imFReader.readLine().split(",");
            fAtBeta[i][1] = Double.parseDouble(imF[1].trim());
        }
        gSeries.rotateFtoG(fAtBeta);
        gSeries.setgAtBeta(fAtBeta);
        //[261.0016582858428, 28.452144305679546, 2.3693797179877887], 261.07309238484356, 
        //[261.31522681873514, -6.883771986248166, 0.08672183144103376]
        //        1.7592847984372848 ** 0.03495629927494104
        final int initialPadding = 40;
        double zeta = evaluateZeta(261.07309238484356, initialPadding, gSeries);
        System.out.println( " zeta " + zeta);
        double[] zero = {261.0016582858428, 261.31522681873514};
        double[] expectedDer = {28.452144305679546, -6.883771986248166};
        for (int i = 0; i < zero.length; i++) {
            validateZero(zero[i], expectedDer[i], initialPadding, gSeries,false);
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
        System.out.println( " zeta " + zeta);
        if(checkAssert){
        assertTrue(Math.abs(zeta) < 0.000001);
        }
        
        double delta = 0.01*gAtBeta.spacing;
        double zetaplus = evaluateZeta(zero+delta, initialPadding, gAtBeta);
        double zetaminus = evaluateZeta(zero-delta, initialPadding, gAtBeta);
        //244.92059950582586, 23.85164367971759, 1.0396728565623998
        double der = (zetaplus-zetaminus)/(2*delta);
        System.out.println("der " + der);
        if(checkAssert){
        assertEquals(expectedDer,der, 0.001);
        }
    }

    public static double evaluateZeta(double zero, final int initialPadding, GSeries gAtBeta) {
        double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, initialPadding, 1.6E-9, false);
        double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
        return zeta;
    }
}
