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
import java.io.PrintStream;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.Arrays;

import org.junit.Ignore;
import org.junit.Test;

import riemann.Gram;
import riemann.Interpolate;
import riemann.Rosser;
import riemann.Interpolate.Poly3;
import riemann.Interpolate.Poly4;
import riemann.PolyInterpolate;
import riemann.Rosser.ZeroInfo;

public class MoreGSeriesTest {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
//        nf.setMinimumIntegerDigits(2);
        nf.setMinimumFractionDigits(4);
        nf.setMaximumFractionDigits(4);
        nf.setGroupingUsed(false);
    }

    @Test  
    public void testSymmetryRelations() throws Exception{
        //check the symmetry and antisymmetry relations from the output of distributions
    	//input can come from zetaHist.R
    	//input for zetaHist.R can come from testInterpolate()
    	int k = 12;
        //File file = new File("data/gzetaE12/real" + k +  ".txt");
        File file = new File("data/gzetaE12/calc12.txt");
        BufferedReader zeroIn = new BufferedReader(new FileReader(file));
        int k2 = 2*k;
        double[][] vals = new double[k2][6];
        for (int i = 0; i < k2; i++) {
            String in = zeroIn.readLine();
            in = in.replace("\\\\", "");
            String[] line = in.split("[&\\s]+");
            for (int j = 0; j < vals[i].length; j++) {
                vals[i][j]= Double.parseDouble(line[j+1].trim());
            }
        }
        zeroIn.close();
        //anti
        System.out.println("antisymmetry");
        System.out.println("q1 median mean q3");
        for (int i = 0; i < k2/2; i++) {
            for (int j = 1; j < vals[i].length-1; j++) {
               double diff = -100;
               if(j == 2 || j ==3){
            	   // median, mean
                   diff = vals[i][j] + vals[i+k2/2][j];
               } else {
                   diff = vals[i][j] + vals[i+k2/2][6-1-j];
               }
                System.out.print(nf.format(diff) + " ");
            }
            System.out.println(" | " + i);
        }
        System.out.println("symmetry");
        System.out.println("q1 median mean q3");
        //symm
        for (int i = 1; i < k2/2; i++) {
            for (int j = 1; j < vals[i].length-1; j++) {
                double diff = -100;
                diff = vals[i][j] - vals[k2-i][j];
                System.out.print(nf.format(diff) + " ");
            }
            System.out.println(" | " + i);
        }
        //median/mode
        System.out.println("median/mode");
        for (int i = 0; i < k2; i++) {
           {
                double diff =  vals[i][2]/ vals[i][3];
                System.out.print(i + " " + 
                        nf.format(diff) + " " + nf.format(vals[i][4]- vals[i][1]));
            }
           System.out.println(" | " + i);
        }
    }
    
    
    @Test @Ignore 
    public void testKurtosisSymmetryRelations() throws Exception{
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

    @Test @Ignore 
    public void testFix() throws Exception{
        GSeries gSeries = Interpolate.readGSeries();
        System.out.println("gSeries.begin " + gSeries.begin);
        final int initialPadding = 40;
        double[][] tmin = {
            {480.82757562193734, -12.479455830100015, 0.13155469450002938, 480.8478824361364},
            {710.5642631336118, 1511.1436494073869, 167.70620377103543, 710.6839237435237,},
            {3811.050548669376, 570.9429277899857, 386.1396790368941, 3811.2260977681094,},
            {90977.64166585186, 644.799901929005, 513.7189446414618, 90977.81218358524,},
            {100798.08697164342, 83.55187028339371,  1.2667769416536376, 100798.11938493361,},
          };
        NumberFormat nf1 = NumberFormat.getInstance();
        nf1.setMinimumFractionDigits(2);
        nf1.setMaximumFractionDigits(2);
        nf1.setGroupingUsed(false);

        for (int i = 0; i < tmin.length; i++) 
        {
           System.out.println(); 
           double d = tmin[i][3];
            double expectedDer = tmin[i][1];
            double expected = expectedDer<0?-tmin[i][2]:tmin[i][2];
            double min = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            System.out.println(nf.format(d) + ", " + nf.format(min) + ", " + nf.format(expected) + ", "
                    + nf1.format((min - expected) * 100.0 / expected) + ", " + nf1.format((min - expected)));

            int midIdx1 = -1;
            double[] g0incr = FixGSeries.evalGSeriesIncrement(gSeries, midIdx1, initialPadding, tmin[i]);
            midIdx1 = gSeries.midIdx;
            gSeries.incrementGValueAtIndex(midIdx1, g0incr);
            System.out.println("eval zero final " + Interpolate.evaluateZeta(tmin[i][0], initialPadding , gSeries) + ", "
                    + Interpolate.evaluateDer(tmin[i][0], initialPadding , gSeries)
                    + " cf " + expectedDer);

            min = Interpolate.evaluateZeta(d, initialPadding, gSeries);
            System.out.println(nf.format(d) + ", " + nf.format(min) + ", " + nf.format(expected) + ", "
                    + nf1.format((min - expected) * 100.0 / expected) + ", " + nf1.format((min - expected)));
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
    	//input can come from riemann.Interpolate.consolidatedF()
    	//That method stores the output G series from riemann.Interpolate.readItems()
    	//provides input for zetaHist.R 
    	// after that is run, use the testSymmetryRelations method
        File gFile = new File("out/gSeries" + Interpolate.prefix + "/gSeriesConsolidated.dat");
        GSeries gSeries = Interpolate.readGSeries(gFile);
        final int initialPadding = 40;
        
//        double[] zero = {100415.50500735927, 100415.61036506912, 
//                100797.8878505715, 100798.08697164342,  };
//        double[] expectedDer = {-46.06567120662985, 45.21334158268663, 
//                -152.8048262150694, 83.55187028339371, };
//        for (int i = 0; i < zero.length; i++) {
//            Interpolate.validateZero(zero[i], expectedDer[i], initialPadding, gSeries,false);
//        }

        int k = 12;
        File file = new File("out/gzeta" + Interpolate.prefix + "/gzeta"
        		+ k + ".csv");
        if (!file.getParentFile().exists()) {
            file.getParentFile().mkdirs();
        }
        PrintWriter out = new PrintWriter(file);
        writeZetaPhi( out, gSeries, initialPadding, 2*gSeries.spacing, k);
        out.close();
    }

    private void writeZetaPhi( PrintWriter out, GSeries gAtBeta, 
            int initialPadding, double incr, int k) throws Exception{
        double[] oddsum = new double[k], evensum = new double[k];
        final double[] zeta = new double[2*k];
        final double firstGram = gAtBeta.begin + initialPadding*incr;
        int N = gAtBeta.gAtBeta.length/2 - 2*initialPadding;
        double gram = firstGram;
        for (int i = 0; i < N; i++) {
            gram += incr;
            for (int j = 0; j < k; j++) {
                double t = gram + j*incr/k;
                double[] gFromBLFI;
                try {
                    gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( 
                            t, 4, initialPadding, 1.6E-9, false);
                } catch (Exception e) {
                    System.out.println("i " + i);
                    throw e;
                }
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
        System.out.println(Arrays.toString(evensum));
        for (int j = 0; j < k; j++) {
            System.out.print(nf.format(oddsum[j]) + " ");
            assertEquals(-2.00*Math.cos(j*Math.PI/k), oddsum[j], 0.1);
        }
        System.out.println();
        for (int j = 0; j < k; j++) {
            System.out.print(nf.format(evensum[j]) + " ");
            assertEquals(2.00*Math.cos(j*Math.PI/k), evensum[j], 0.1);
        }
        System.out.println();
    }
    
    @Test @Ignore 
    public void testGenerateE12() throws Exception{
    	double[][] zetaGramMean = {{0, 0},{0, 0}};
    	double[][] zetaMidMean = {{0, 0},{0, 0}};
    	int noffset = Interpolate.noffset;
        int correction = -1;
		double baseLimit = Interpolate.baseLimit;
		double begin= baseLimit  + 
        		(noffset-correction )* (Interpolate.gramIncr);
        double gramIncr = Interpolate.gramIncr;
		GSeries gSeries = new GSeries(1, 0, Interpolate.offset, begin, gramIncr );
        
        double zetaCorrection1 = GSeries.correction( gSeries.basesqrtArg1);
        PrintStream out = Interpolate.out;
		BufferedReader[] zeroIn = Interpolate.zeroIn;
		Interpolate.zeroInput = Rosser.readZeros(baseLimit, out , zeroIn , null);
        System.out.println(Arrays.toString(Interpolate.zeroInput.lastZero)  +
                ", " + baseLimit + ",\n " + Arrays.toString(Interpolate.zeroInput.nextValues));
        
        String zetaFile = "data/zetaE12.csv";
        BufferedReader zetaIn = 
                new BufferedReader(new FileReader(zetaFile));
        zetaIn.readLine();
        for (int i = 0; i < 2-correction; i++) {
            zetaIn.readLine();
		}

        long sampleSize = 1000100;
        int count = 0;
        Poly3 poly = null;
		while (count  < sampleSize  ) {
            int n = count + noffset;
            double upperLimit = baseLimit + (n-correction-1)* (gramIncr);
            poly = updateZero(zetaGramMean, zetaCorrection1, 
            		out, zeroIn, poly, n, upperLimit);
            upperLimit += gramIncr/2;
            poly = updateZero(zetaMidMean, zetaCorrection1, 
            		out, zeroIn, poly, n, upperLimit);
            count++;
        }    
        zetaIn.close();
        System.out.println("*** zetaGram_MeanOdd poly " + nf.format(2*zetaGramMean[1][0]/sampleSize));
        System.out.println("*** zetaGram_MeanEven poly " + nf.format(2*zetaGramMean[0][0]/sampleSize));
        System.out.println("*** zetaGram_MeanOdd poly5 " + nf.format(2*zetaGramMean[1][1]/sampleSize));
        System.out.println("*** zetaGram_MeanEven poly5 " + nf.format(2*zetaGramMean[0][1]/sampleSize));

        System.out.println("*** zetaMid_MeanOdd poly " + nf.format(2*zetaMidMean[1][0]/sampleSize));
        System.out.println("*** zetaMid_MeanEven poly " + nf.format(2*zetaMidMean[0][0]/sampleSize));
        System.out.println("*** zetaMid_MeanOdd poly5 " + nf.format(2*zetaMidMean[1][1]/sampleSize));
        System.out.println("*** zetaMid_MeanEven poly5 " + nf.format(2*zetaMidMean[0][1]/sampleSize));
}


	private Poly3 updateZero(double[][] zetaGramMean, double zetaCorrection1, PrintStream out, BufferedReader[] zeroIn,
			Poly3 poly, int n, double upperLimit) throws FileNotFoundException, IOException {
		double z0;
		double z1;
		double d0;
		double d1;
		double max;
		if(upperLimit<=Interpolate.zeroInput.nextValues[0]){
			Interpolate.zeroInput = new ZeroInfo(0, Interpolate.zeroInput);
		} else {
			Interpolate.zeroInput = Rosser.readZeros(upperLimit , out, zeroIn,  
					Interpolate.zeroInput.nextValues);
		    z0 = Interpolate.zeroInput.lastZero[0];
		    z1 = Interpolate.zeroInput.nextValues[0];
		    d0 = Interpolate.zeroInput.lastZero[1];
		    d1 = Interpolate.zeroInput.nextValues[1];
		    max = d0>0?Interpolate.zeroInput.lastZero[2]:-Interpolate.zeroInput.lastZero[2];
		    Interpolate.poly = new Poly4(z0,z1, d0,d1,max);
		    
		    ZeroPoly zeroPoly = new ZeroPoly(Rosser.zeros, Rosser.derivatives);
		    double secondDer = zeroPoly.secondDer(1);
		    poly = new PolyInterpolate.Poly5(z0,z1, d0,d1,secondDer,max);
		}
		double zeta = Interpolate.poly.eval(upperLimit)- zetaCorrection1;
//            System.out.println("gram zeta " + zeta + ", " +  upperLimit + " n upperLimit (" + (n+1) +")");
//            String input = zetaIn.readLine();
//            String[] parsed = input.split(",");
//            System.out.println(Arrays.toString(parsed));
		final int nmod2 = n%2;
		zetaGramMean[nmod2][0] += zeta;
		
		double zeta5 = poly.eval(upperLimit)- zetaCorrection1;
		zetaGramMean[nmod2][1] += zeta5;
		return poly;
	}
    
    @Test @Ignore 
    public void testEstimateE12() throws Exception{
    	double[] zetaGramMean = new double[]{0, 0};
        int count = 0;
        int N = 1000082;
        String zetaFile = "data/zetaE12.csv";
        BufferedReader zetaIn = 
                new BufferedReader(new FileReader(zetaFile));
        zetaIn.readLine();
        double base = Interpolate.baseLimit;
        for (int i = 0; i < 41; i++) {
            zetaIn.readLine();
		}

        GSeries gSeries1 = Interpolate.readGSeries();
//        long sampleSize = N-2*Interpolate.initialPadding;
        long sampleSize = 1005;
        int deviations = 0;
        while (count < sampleSize  ) {
            int n = count + Interpolate.noffset+Interpolate.initialPadding;
            double t = base + (n-2)* (Interpolate.gramIncr);
            double zeta =  Interpolate.evaluateZeta(t, Interpolate.initialPadding, gSeries1);        
            String input = zetaIn.readLine();
            String[] parsed = input.split(",");
            double zetaSaved =  Double.parseDouble(parsed[1]);  
            double diff = Math.abs(zeta-zetaSaved);
            if(diff > 0.1) {
            	deviations++;
            	System.out.println( n + " " + nf.format(t) + " " + nf.format(zeta)
            		+ " " + nf.format(diff));
            }
            final int nmod2 = n%2;
            zetaGramMean[nmod2] += zeta;
            count++;
        }  
        zetaIn.close();
        System.out.println("deviations " + deviations);
//        System.out.println("*** zetaGram_MeanOdd " + 2*zetaGramMean[1]/sampleSize);
//        System.out.println("*** zetaGram_MeanEven " + 2*zetaGramMean[0]/sampleSize);
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
