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
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.function.Function;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import math.GSeries;
import math.ZeroPoly;
import riemann.Rosser.ZeroInfo;

public class Interpolate {
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
    static int breaks = 0;
    static double zetaCorrection1;
    static  double absMax = 0;
    public static int correction = 0;
    public static double baseLimit;
    public static double gramIncr;
    public static int noffset;
    public static String prefix;
    public static final int initialPadding = 40;

    
    public abstract static class Poly3{
        final double z0, z1;
        final double d0, d1;
        double denom;
        final double h;
        
        final double pow2; 
        final double pow1;
        final double pow0;
        final double cross;
        
        Poly3(double z0, double z1, double d0, double d1) {
            this.z0 = z0;
            this.z1 = z1;
            this.d0 = d0;
            this.d1 = d1;
            h = (z1-z0);
            denom = h*h;
            
            pow2 = (d1+d0)/denom;
            pow0 = ((d1+d0)*(2*z0*z1) + z0*z0*d1  +z1*z1*d0)/denom;
            cross = (z0*d1+z1*d0)/denom;
            pow1 = 2*(pow2*(-z0-z1) - cross);
        }
        public double eval1(double x){
            double ret = (x-z0)*(x-z1)/denom;
            ret *= ((x-z0)*d1 + (x-z1)*d0);
//            double ret = (x-z0)*(x-z1);
//            ret *= pow2*x - cross;
            return ret;
        }
        
        public double der(double x){
            double ret = (
                  (d1+d0)*(x-z0)*(x-z1) + (x-z0)*((x-z0)*d1 
                  + (x-z1)*d0) + (x-z1)*((x-z0)*d1 + (x-z1)*d0)
                               )/denom;
//            double ret = (3*x*x*pow2 
//                    + x*pow1
//                    + pow0);
            return ret;
        }

        public double secondDer(double x){
            double ret = (6*x*pow2 
                    + pow1);
            return ret;
        }

        public double secondDerRHS(){
            double ret = 2*(d0+2*d1)/h;
            return ret;
        }
        
        abstract void estimateC(  double xmin) ;
        public abstract double eval(double x);
        public abstract double getPositionMax();
       
        protected double processMax( ) {
            double[] oldest = new double[]{z0, d0, 0};
            double[] upper = new double[]{z1, d1, 1};
            double[] wts = new double[]{Math.abs(d1),Math.abs(d0)};
            double xmin = z0 - d0*(z1-z0)/(d1-d0);
            double precision = 0.001*Math.abs(d0);
            for (int i = 0; i < 10; i++) {
               Interpolate.xmin(oldest, upper, wts, precision,  (x)->der(x));
               if(oldest[2]> 99){
                   xmin = oldest[0];
                   break;
               }
               wts[0] = Math.abs(upper[1]);
               wts[1] = Math.abs(oldest[1]);
               if(wts[0]>wts[1]){
                   xmin = oldest[0];
               }else {
                   xmin = upper[0];
               }
            }
            estimateC(xmin);
            oldest = new double[]{z0, d0, 0};
            upper = new double[]{z1, d1, 1};
            double dermin = der(xmin);
            Interpolate.fixLimits(oldest, upper, xmin, dermin);            
            wts = new double[]{Math.abs(upper[1]),Math.abs(oldest[1])};
            precision = 0.0000001;
            for (int i = 0; i < 10; i++) {
               estimateC(xmin);
               Interpolate.xmin(oldest, upper, wts, precision, (x)->der(x));
               if(oldest[2]> 99){
                   xmin = oldest[0];
                   break;
               }
               wts[0] = Math.abs(upper[1]);
               wts[1] = Math.abs(oldest[1]);
               if(wts[0]>wts[1]){
                   xmin = oldest[0];
               }else {
                   xmin = upper[0];
               }
            }
            return xmin;
        }

    };

    public static class Poly4 extends Poly3{
        double C;
        double positionMax;

        @Override
        public double getPositionMax() {
			return positionMax;
		}
		double max;
        double cdenom;
        final double sum2;
        public Poly4(double z0, double z1, double d0, double d1, double max) {
            super(z0, z1, d0, d1);
            this.max= max;
            positionMax = processMax();
            cdenom = 4*C/denom;
            sum2 = (z1+z0)/2.0;
        }
        
        void estimateC(  double xmin) {
            double mult = (xmin-z0)*(xmin-z1);
            mult = mult*mult/denom;
            C = (max - eval1(xmin))/mult;
            cdenom = 4*C/denom;
        } 

        public double eval(double x){
            double mult = (x-z0)*(x-z1);
            mult = mult*mult/denom;
            double ret = super.eval1(x)+C*mult;
            return ret;
        }        
        public double der(double x){
            double ret = super.der(x);
            //ret += (x-z0)*(x-z1)*(x-sum2)*cdenom;
            ret += 2*C*(x-z0)*(x-z1)*(2*x-z1-z0)/denom;
            return ret;
        }
        public double secondDer(double x){
            double ret = super.secondDer(x);
            ret += ((x-z0)*(x-z1) + 2*(x-sum2)*(x-sum2))*cdenom;
            return ret;
        }
        public double secondDerRHS(){
            double ret = super.secondDerRHS();
            return ret += 2*C;        
        }
    }


    static void xmin(double[] oldest, double[] upper, 
            double[] wts, double precision, 
            Function<Double, Double> derivativeFunction){
        double x0 = oldest[0];
        double x1 = upper[0];
        double xmin = (wts[0]*x0 + wts[1]*x1)/(wts[0]+wts[1]);
        double dermin = derivativeFunction.apply(xmin);
        if(Math.abs(dermin)<precision){
            oldest[0]=xmin;
            oldest[1]=dermin;
            oldest[2]=100;
            return;
        } 
        Interpolate.fixLimits(oldest, upper, xmin, dermin);
        xmin = (oldest[0] + upper[0])/(2);
        dermin = derivativeFunction.apply(xmin);
        Interpolate.fixLimits(oldest, upper, xmin, dermin);
    }

    static void fixLimits(double[] oldest, double[] upper, double xmin, double dermin) {
        if(Math.signum(dermin) == Math.signum(oldest[1])){
            oldest[0]=xmin;
            oldest[1]=dermin;
        } else {
            upper[0]=xmin;
            upper[1]=dermin;
        }
    }
    
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
            double zeta =  getZetaEstimate(nprime, idx, upperLimit, zetaGramMean, fAtBeta, 0);

//            System.out.println();
//            System.out.println(Arrays.toString(zeroInput.lastZero)  +
//                    ", \n" + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
//            System.out.println("gram 2*zeta " + 2*zeta + ", " +  upperLimit + " n upperLimit (" + (n+1) +")");
            
            upperLimit += gramIncr/2;
            zeta = getZetaEstimate(nprime, idx, upperLimit, zetaMidMean,imFmid,1);
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
        double[][] consolidated = consolidatedF();
        ;
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

    private static double getZetaEstimate(int nprime, int idx, double upperLimit, 
            double[] zetaMean, double[][] fStorageReIm, int gramOrMid) throws Exception {
    	double zetaEstMid = updateZeroInput(upperLimit);
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
        //i == 0, Gram
        switch(nmod2){
        case 0:
            fStorageReIm[idx][gramOrMid] = (-zeta);
            if(gramOrMid==0){
                //i == 0, Gram
                gramDer[idx][0] = -poly.der(upperLimit);
                gramDer[idx][1] = -poly.secondDer(upperLimit);
            }
            break;
        case 1:
            fStorageReIm[idx][gramOrMid] = (zeta);
            if(gramOrMid==0){
                gramDer[idx][0] = poly.der(upperLimit);
                gramDer[idx][1] = poly.secondDer(upperLimit);
            }
        }
        if(idx>0 && gramOrMid==0){
            imFmid[idx-1][0] = -100;
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
    private static final double updateZeroInput(double upperLimit) throws FileNotFoundException, IOException {
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
//            if(Double.isFinite(Rosser.zeros[0])) {
//	            ZeroPoly zeroPoly = new ZeroPoly(Rosser.zeros, Rosser.derivatives);
//	            double secondDer = zeroPoly.secondDer(1);
//	            poly = new PolyInterpolate.Poly5(z0,z1, d0,d1,secondDer,max);
//            	
//            } else {
            	poly = new Poly4(z0,z1, d0,d1,max);
//            }
        }
        double zetaEstMid = poly.eval(upperLimit);
        return zetaEstMid;
    }
    
    public static double[][]  consolidatedF(  ) throws IOException {
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
        
        double[][] consolidated = new double[2*imFmid.length][2];
        return consolidated;
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
