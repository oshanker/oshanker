package riemann;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import math.GSeries;
import math.ZeroPoly;
import riemann.Rosser.ZeroInfo;

public class Interpolate {
	public enum PolyOption{
		USE_POLY4,
		USE_MIXED,
		USE_POLY5,
	};
	
	public enum GramOrMid{
		GRAM,
		MID
	}
	
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
        getZerosFile();
    }
    public static BigDecimal offset;
    public final static PrintStream out = null;
    static PrintStream validateOut = null;
    public static BufferedReader[] zeroIn;
    static double[] lastZeroSeen1;
    static double[][] gramDer;
    static double[][] fAtBeta;
    static double[][] imFmid;
    
    public static ZeroInfo zeroInput;
    public static Poly3 poly = null;
    public static Poly3 poly5 = null;
    public static PolyOption polyOption = PolyOption.USE_MIXED;
    
    static int breaks = 0;
    static double zetaCorrection1;
    static  double absMax = 0;
    public static int correction = 0;
    public static double baseLimit;
    public static double gramIncr;
    public static int noffset;
    public static String prefix;
    public static final int initialPadding = 40;


    static void getZerosFile()  {
        try {
            Rosser.readConfig("data/RosserConfig.txt");
            String zerosFile = Rosser.getParam("zerosFile");
            System.out.println("zerosFile " + zerosFile);
            zeroIn = new BufferedReader[] { 
                    new BufferedReader(new FileReader(zerosFile)), 
                    null,
                    null 
                    };

            String derFile = zerosFile + ".der";
            zeroIn[1] = new BufferedReader(new FileReader(derFile));
            String maxFile = zerosFile + ".max";
            zeroIn[2] = new BufferedReader(new FileReader(maxFile));
            lastZeroSeen1 = new double[zeroIn.length + 1];
            offset = new BigDecimal(Rosser.getParam("bdoffset"));
            correction = 0;
            if(Rosser.configParams.containsKey("correction")){
                correction = Rosser.getParamInt("correction");
            }
            baseLimit = Rosser.getParamDouble("baseLimit");
            gramIncr = Rosser.getParamDouble("gramIncr");
            noffset = Rosser.getParamInt("noffset");
            Pattern p = Pattern.compile("\\S*stats(E\\d\\d)\\S*");
            Matcher m = p.matcher(Rosser.getParam("conjecturesOutFile"));
            m.matches();
            prefix = m.group(1);
            
        } catch (Exception e){
            throw new IllegalStateException(e);
        }
    }

    private static void readItems(   )
            throws Exception {
        double begin= baseLimit + (noffset-correction)* (gramIncr);
        GSeries gSeries1 = new GSeries(1, 0, offset, begin, gramIncr);
        
        zetaCorrection1 = GSeries.correction( gSeries1.basesqrtArg1);
        
        BigDecimal tvalsi = offset.add(BigDecimal.valueOf(begin), Gram.mc);
        BigDecimal gramIndex1 = Gram.theta(tvalsi, Gram.mc).divide(Gram.pi, Gram.mc);
        String[] line = Rosser.getParam("header").split("[-L]+");
        
        gramIndex1 = 
                gramIndex1.subtract(new BigDecimal(line[1]), Gram.mc);
        System.out.println( gSeries1.begin + ", zetaCorrection " + zetaCorrection1
                + ", gram index " + gramIndex1 + ", " + line[1]);
        
        int N = Rosser.getParamInt("N");
        //        N = 12;
        N = 1000102;
        int count = 0;
        zeroInput = Rosser.readZeros(baseLimit, out, zeroIn, null);
        System.out.println(Arrays.toString(zeroInput.lastZero)  +
                ", " + baseLimit + ", " + Arrays.toString(zeroInput.nextValues));
        System.arraycopy(zeroInput.nextValues, 0, lastZeroSeen1, 0, zeroIn.length);
        //f at Mid. 0 -> real. 1 -> im
        imFmid = new double[N][2];
        //f at Gram. 0 -> real. 1 -> im
        fAtBeta = new double[N][2];
        gramDer = new double[N][2];
        double[] zetaMidMean = {0, 0};
        double[] zetaGramMean = {0, 0};
        int idx = 0;
        while (count < N  ) {
        	//this n is actually n-1!!!
        	//idx = 0, n = 3 in zetaE12.csv
            int nprime = count + noffset;
            double upperLimit = baseLimit + (nprime-correction-1)* (gramIncr);
            // populate fAtBeta,  zetaGramMean
            double zeta =  getZetaEstimate(nprime, idx, upperLimit, zetaGramMean, fAtBeta,  GramOrMid.GRAM);

//            System.out.println();
//            System.out.println(Arrays.toString(zeroInput.lastZero)  +
//                    ", \n" + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
//            System.out.println("gram 2*zeta " + 2*zeta + ", " +  upperLimit + " n upperLimit (" + (n+1) +")");
            
            upperLimit += gramIncr/2;
            zeta = getZetaEstimate(nprime, idx, upperLimit, zetaMidMean,imFmid, GramOrMid.MID);
            //            System.out.println("mid  2*zeta " + 2*zeta + ", " +  upperLimit + " (" + (n+1) +")");
            if (count == N - 1) {
                System.out.println("final n " + nprime);
                System.out.println();
                System.out.println(Arrays.toString(zeroInput.lastZero) + ", \n" + upperLimit + ", "
                        + Arrays.toString(zeroInput.nextValues));
                System.out.println("mid " + zeta + ", " +  upperLimit + " (" + (nprime+1) +")");
            }
            idx++;
            count++;
        }
        System.out.println( "breaks: " + breaks);
        System.out.println("*** zetaMidMeanOdd " + zetaMidMean[1]/N);
        System.out.println("*** zetaGramMeanOdd " + zetaGramMean[1]/N);
        System.out.println("*** zetaMidMeanEven " + zetaMidMean[0]/N);
        System.out.println("*** zetaGramMeanEven " + zetaGramMean[0]/N);
        
//        double predictedSqrtArg1 = gSeries.basesqrtArg1 + gSeries.dsqrtArg1*N*gramIncr;
//        zetaCorrection1 = GSeries.correction( predictedSqrtArg1);
//        System.out.println( "final zetaCorrection: " + zetaCorrection1);
        //imFGramPoints( );
        //reFMidGramPoints();
        consolidatedF();
        double[][] consolidated = new double[2*imFmid.length][2];
		for (int i = 0; i < imFmid.length; i++) {
            consolidated[2*i][0] = fAtBeta[i][0];
            consolidated[2*i][1] = fAtBeta[i][1];
            consolidated[2*i+1][0] = imFmid[i][0];
            consolidated[2*i+1][1] = imFmid[i][1];
        }
        begin= baseLimit + (noffset-correction-1)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr/2);
        System.out.println( "gSeries begin: " + gSeries.begin + ", " );

        gSeries.rotateFtoG(consolidated);
        consolidated[0][1] =  Double.NEGATIVE_INFINITY;  
        gSeries.setgAtBeta(consolidated);
        File gFile = new File("out/gSeries" + prefix + "/gSeriesConsolidated.dat");
        storeG(gSeries.begin, gSeries.spacing, gSeries.gAtBeta, gFile);
        /*
        # python -i plot_distribution.py 
        # load dataset
        # input generated in math.GSeriesTest.testWriteZetaPhiE12()
        # or in math.MoreGSeriesTest.testInterpolate()
        */
        zetaMidMean = new double[]{0, 0};
        zetaGramMean = new double[]{0, 0};
        count = 0;
        GSeries gSeries11 = readGSeries();
        double base = baseLimit + gramIncr/4;
        long sampleSize = N-2*initialPadding;
        while (count < sampleSize  ) {
            int n = count + noffset+initialPadding;
            double upperLimit = base + (n-correction-1)* (gramIncr);
            double zeta =  Interpolate.evaluateZeta(upperLimit, initialPadding, gSeries11);        
            final int nmod2 = n%2;
            zetaGramMean[nmod2] += zeta;
            count++;
        }   
        System.out.println("*** zetaGram_MeanOdd " + 2*zetaGramMean[1]/sampleSize);
        System.out.println("*** zetaGram_MeanEven " + 2*zetaGramMean[0]/sampleSize);
        
//        for (int i = 0; i < 2; i++) {
//            for (int j = 0; j < 2; j++) {
//                System.out.println("fAtBeta[" +  i + "," + j + "] = " + fAtBeta[i][j]);
//                System.out.println("imFmid[" +  i + "," + j + "] = " + imFmid[i][j]);
//            }
//        }
//        for (int i = fAtBeta.length-2; i < fAtBeta.length; i++) {
//            for (int j = 0; j < 2; j++) {
//                System.out.println("fAtBeta[" +  i + "," + j + "] = " + fAtBeta[i][j]);
//                System.out.println("imFmid[" +  i + "," + j + "] = " + imFmid[i][j]);
//            }
//        }
    }

    public static double getZetaEstimate(int nprime, int idx, double upperLimit, 
            double[] zetaMean, double[][] fStorageReIm, GramOrMid gramOrMidEnum) throws Exception {
    	double zetaEstMid = updateZeroInput(upperLimit, gramOrMidEnum);
        if(Math.abs(zeroInput.lastZero[2])>absMax){
            absMax = Math.abs(zeroInput.lastZero[2]);
            if(absMax>130){
                System.out.println();
                System.out.println(Arrays.toString(zeroInput.lastZero)  +
                        ", \n positionMax " + (poly instanceof Poly4?poly.getPositionMax() :"unknown") +
                        ", der " + poly.der(poly.getPositionMax()) + 
                        ", \n" + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
                System.out.println("secondDerRHS " + poly.secondDerRHS() 
                + ", zetaEstMid " + zetaEstMid + " (" + (nprime+1) +")");
             }
        }
        double zeta = (zetaEstMid - zetaCorrection1)/2;
        final int nmod2 = nprime%2;
        // mean of zeta
        zetaMean[nmod2] += zeta;
        int gramOrMidIndex;
        switch (gramOrMidEnum) {
		case GRAM:
			gramOrMidIndex = 0;
			break;

		default:
			gramOrMidIndex = 1;
			break;
		}
        
		//i == 0, Gram
        switch(nmod2){
        case 0:
            fStorageReIm[idx][gramOrMidIndex] = (-zeta);
            if(gramOrMidIndex==0){
                //i == 0, Gram
                gramDer[idx][0] = -poly.der(upperLimit);
                gramDer[idx][1] = -poly.secondDer(upperLimit);
            }
            break;
        case 1:
            fStorageReIm[idx][gramOrMidIndex] = (zeta);
            if(gramOrMidIndex==0){
                gramDer[idx][0] = poly.der(upperLimit);
                gramDer[idx][1] = poly.secondDer(upperLimit);
            }
        }
        return zeta;
    }

    /**
     * should we pass Gram/Mid info?
     * @param upperLimit
     * @return 
     * @throws FileNotFoundException
     * @throws IOException
     */
    private static final double updateZeroInput(double upperLimit, GramOrMid gramOrMid
    		) throws FileNotFoundException, IOException {
        if(upperLimit<=zeroInput.nextValues[0]){
            zeroInput = new ZeroInfo(0, zeroInput);
        } else {
            zeroInput = Rosser.readZeros(upperLimit , out, zeroIn,  
                    zeroInput.nextValues);
            if(lastZeroSeen1[0] != zeroInput.lastZero[0]){
                breaks++;
            } 
            System.arraycopy(zeroInput.nextValues, 0, lastZeroSeen1, 0, zeroIn.length);
            final double z0 = zeroInput.lastZero[0];
            final double z1 = zeroInput.nextValues[0];
            final double d0 = zeroInput.lastZero[1];
            final double d1 = zeroInput.nextValues[1];
            final double max = d0>0?zeroInput.lastZero[2]:-zeroInput.lastZero[2];
            switch (polyOption) {
			case USE_POLY4:
				poly = new Poly4(z0,z1, d0,d1,max);
				break;

			case USE_MIXED:
				poly = new Poly4(z0,z1, d0,d1,max);
	            if(Double.isFinite(Rosser.zeros[0])) {
		            ZeroPoly zeroPoly = new ZeroPoly(Rosser.zeros, Rosser.derivatives);
		            double secondDer = zeroPoly.secondDer(1);
		            poly5 = new PolyInterpolate.Poly5(z0,z1, d0,d1,secondDer,max);
	            	
	            } else {
	            	poly5 = new Poly4(z0,z1, d0,d1,max);
	            }
				break;

			default:
				//use poly5
	            if(Double.isFinite(Rosser.zeros[0])) {
		            ZeroPoly zeroPoly = new ZeroPoly(Rosser.zeros, Rosser.derivatives);
		            double secondDer = zeroPoly.secondDer(1);
		            poly = new PolyInterpolate.Poly5(z0,z1, d0,d1,secondDer,max);
	            	
	            } else {
	            	poly = new Poly4(z0,z1, d0,d1,max);
	            }
				break;
			}
        }
        double zetaEstMid;
        switch (polyOption) {
		case USE_MIXED:
	        switch (gramOrMid) {
			case GRAM:
    			zetaEstMid = poly5.eval(upperLimit);
				break;

			default:
    			zetaEstMid = poly.eval(upperLimit);
				break;
			}
			break;
		default:
			zetaEstMid = poly.eval(upperLimit);
			break;
        
        }
       return zetaEstMid;
    }
    
    public static void  consolidatedF(  ) throws IOException {
//        File file = new File(Rosser.getParam("conjecturesOutFile")
//                .replace("stats", "validateConsolidatedF"));
//        validateOut = new PrintStream(file);
        double[] yF = new double[imFmid.length];
        for (int i = 0; i < yF.length; i++) {
            yF[i] = imFmid[i][1];
        }
        
        NormalizedSpline normalizedSplineF = new NormalizedSpline(yF);
        normalizedSplineF.evalMid(fAtBeta, 1, 1);

        double[] yIm = new double[imFmid.length];
        for (int i = 0; i < yIm.length; i++) {
            yIm[i] = fAtBeta[i][0];
        }
        
        NormalizedSpline normalizedSplineIm = new NormalizedSpline(yIm);
        normalizedSplineIm.evalMid(imFmid, 0, 0);
        
//        readAndValidate();
//        validateOut.close();
    }
        
    private static void imFGramPoints(  ) throws IOException {
        File file = new File(Rosser.getParam("conjecturesOutFile")
                .replace("stats", "validate"));
        validateOut = new PrintStream(file);
        double[] y = new double[imFmid.length];
        for (int i = 0; i < y.length; i++) {
            y[i] = imFmid[i][1];
        }
        
        int seriesOffset = 1, position = 1;
        NormalizedSpline normalizedSpline = new NormalizedSpline(y);
        normalizedSpline.evalMid(fAtBeta, seriesOffset, position);

        double begin= baseLimit + (noffset-correction-1)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        System.out.println( "begin: " + gSeries.begin + ", " );

        gSeries.rotateFtoG(fAtBeta);
        fAtBeta[0][1] =  Double.NEGATIVE_INFINITY;  
        gSeries.setgAtBeta(fAtBeta);
        checkZeros(gSeries);
        checkMax(gSeries);
        validateOut.close();
        //        storeG(begin, gramIncr, fAtBeta);
//        readAndValidate();
    }
    
    private static void reFMidGramPoints(  ) throws IOException {
        File file = new File(Rosser.getParam("conjecturesOutFile")
                .replace("stats", "validateMid"));
        validateOut = new PrintStream(file);
        double[] yIm = new double[imFmid.length];
        for (int i = 0; i < yIm.length; i++) {
            yIm[i] = fAtBeta[i][0];
        }
        
        int seriesOffset = 0, position = 0;
        NormalizedSpline normalizedSplineIm = new NormalizedSpline(yIm);
        normalizedSplineIm.evalMid(imFmid, seriesOffset, position);

        double begin= baseLimit + (noffset-correction-1+0.5)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        System.out.println( "begin: " + gSeries.begin + ", " );

        gSeries.rotateFtoG(imFmid);
        gSeries.setgAtBeta(imFmid);
        checkZeros(gSeries);
        checkMax(gSeries);
        validateOut.close();
        //        storeG(begin, gramIncr, fAtBeta);
//        readAndValidate();
    }

    public static GSeries readGSeries() throws FileNotFoundException, IOException {

        File file = new File("out/gSeries" + prefix + "/gSeriesConsolidated.dat");
        return readGSeries(file);
    }

    public static GSeries readGSeries(File file) throws FileNotFoundException, IOException {
        DataInputStream in = dataInputStream(file);
        double begin = in.readDouble();
        double gincr = in.readDouble();
        int R = in.readInt();
        double[][] gBeta = new double[R][2];
        for (int i = 0; i < gBeta.length; i++) {
            gBeta[i][0] = in.readDouble();
            gBeta[i][1] = in.readDouble();
        }
        GSeries gAtBeta = new GSeries(1, 0, offset,  begin,  gincr, gBeta);
        in.close();
        return gAtBeta;
    }

    private static void readAndValidate() throws FileNotFoundException, IOException {
        GSeries gSeries = readGSeries();
        checkZeros(gSeries);
    }

    private static void checkZeros(GSeries gSeries) {
        double[] zero1 = {480.82757562193734, 
                90977.64166585186, 90977.97641173516, 
                100415.50500735927, 100415.61036506912, 
                100797.8878505715, 100798.08697164342,  };
        double[] expectedDer1 = {-12.479455830100015, 
                644.799901929005, -1487.416968799948,
                -46.06567120662985, 45.21334158268663, 
                -152.8048262150694, 83.55187028339371, };
        for (int i = 0; i < zero1.length; i++) {
            validateZero(zero1[i], expectedDer1[i], initialPadding, gSeries,false);
        }
    }

    private static void storeG(double begin, double incr, double[][] gAtBeta) throws IOException, FileNotFoundException {
        File file = new File("out/gSeries" + prefix + "/gSeries.dat");
        storeG(begin, incr, gAtBeta, file);
    }

    private static void storeG(double begin, double incr, double[][] gAtBeta, File file)
            throws FileNotFoundException, IOException {
        DataOutputStream out = outputStream( file);
        out.writeDouble(begin);
        out.writeDouble(incr);
        out.writeInt(gAtBeta.length);
        for (int i = 0; i < gAtBeta.length; i++) {
            out.writeDouble(gAtBeta[i][0]);
            out.writeDouble(gAtBeta[i][1]);
        }
        out.close();
    }

    public static DataInputStream dataInputStream(File file) throws FileNotFoundException {
        InputStream is = new FileInputStream(file);
        // create buffered input stream.
        BufferedInputStream bis = new BufferedInputStream(is);
        // create data input stream to read data in form of primitives.
        DataInputStream in = new DataInputStream(bis);
        return in;
    }

    public static DataOutputStream outputStream(File file) throws FileNotFoundException {
        DataOutputStream out;
        if (!file.getParentFile().exists()) {
            file.getParentFile().mkdirs();
        }
            OutputStream os = new FileOutputStream(file);
            BufferedOutputStream bos = new BufferedOutputStream(os);
            // create data output stream
            out = new DataOutputStream(bos);
        return out;
    }

    public static void validateMax(double postionMax, double expectedMax, 
            final int initialPadding, 
            GSeries gAtBeta, boolean checkAssert) {
        double zeta = Interpolate.evaluateZeta(postionMax, initialPadding, gAtBeta);
        System.out.println( "(max) zeta " + nf.format(zeta) + " cf " + expectedMax
                + " postionMax " + postionMax );
        if(checkAssert){
           if(Math.abs(zeta) > 0.000001){
               throw new IllegalStateException("Expected " + expectedMax + " but got "
                       + zeta);
           }
        }
        double der = evaluateDer(postionMax, initialPadding, gAtBeta);
        System.out.println("der " + nf.format(der) + " cf 0" );
        if(validateOut != null){
            validateOut.println(" postionMax, " + postionMax
                    + ", " + nf.format(zeta)
                    + ", <-, " + nf.format(expectedMax)
                    + ", " + nf.format(der) 
                    + ", " + nf.format(Math.abs(expectedMax-zeta)) 
                    );
        }
        if(checkAssert){
            if(Math.abs(der)>0.001){
                throw new IllegalStateException("Expected 0 "  + " but got "
                        + der);
            }
        }
    }
    
    public static void validateZero(double zero, double expectedDer, 
            final int initialPadding, 
            GSeries gAtBeta, boolean checkAssert) {
        System.out.println( "# zero " + zero );
        double zeta = Interpolate.evaluateZeta(zero, initialPadding, gAtBeta);
        System.out.println( "zeta " + nf.format(zeta) + " cf 0.0" );
        if(checkAssert){
           if(Math.abs(zeta) > 0.000001){
               throw new IllegalStateException("Expected zero but got "
                       + zeta);
           }
        }
        
        double der = evaluateDer(zero, initialPadding, gAtBeta);
        System.out.println("der " + nf.format(der) + " cf " 
        + nf.format(expectedDer));
        if(validateOut != null){
            validateOut.println(" zero, " + zero
                    + ", " + nf.format(zeta)
                    + ", ->, " + nf.format(expectedDer)
                    + ", " + nf.format(der) 
                    + ", " + nf.format(Math.abs(expectedDer-der)) 
                    );
        }
        if(checkAssert){
            if(Math.abs(expectedDer-der)>0.001){
                throw new IllegalStateException("Expected " + expectedDer + " but got "
                        + der);
            }
        }
    }

    public static double evaluateDer(double zero, final int initialPadding, GSeries gAtBeta) {
        double delta = 0.001*gAtBeta.spacing;
        double zetaplus = Interpolate.evaluateZeta(zero+delta, initialPadding, gAtBeta);
        double zetaminus = Interpolate.evaluateZeta(zero-delta, initialPadding, gAtBeta);
        double der = (zetaplus-zetaminus)/(2*delta);
        return der;
    }

    public static double evaluateZeta(double zero, final int initialPadding, GSeries gAtBeta) {
        double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, 
                initialPadding, 1.6E-9, false);
        double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
        return zeta;
    }

    public static void main(String[] args) throws Exception{
        //checkMax();
        readItems();
    }

    private static void checkMax() throws FileNotFoundException, IOException {
        GSeries gSeries = readGSeries();
        checkMax(gSeries);
    }

    private static void checkMax(GSeries gSeries) {
        double[] zero1 = { 682.7988048955597, 3811.2260977681094, 
                66093.42494812438, 90977.81218358524, 
                };
        double[] expectedDer1 = { 138.61973697485362, 386.1396790368941, 
                392.08238609990934, 513.7189446414618
                };
        for (int i = 0; i < zero1.length; i++) {
            validateMax(zero1[i], expectedDer1[i], initialPadding , 
                    gSeries, false);
        }
    }
}
