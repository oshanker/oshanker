package riemann;

import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.Arrays;

import math.GSeries;
import riemann.Rosser.ZeroInfo;
import riemann.PolyInterpolate.Poly5;

public class Interpolate {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
        getZerosFile();
    }
    final static PrintStream out = null;
    static BufferedReader[] zeroIn;
    static double[] lastZeroSeen1;
    static ZeroInfo zeroInput;
    static Poly3 poly = null;
    static int breaks = 0;
    static double zetaCorrection1;
    static         double absMax = 0;

    
    public abstract static class Poly3{
        final double z0, z1;
        final double d0, d1;
        double denom;
        final double h;
        Poly3(double z0, double z1, double d0, double d1) {
            this.z0 = z0;
            this.z1 = z1;
            this.d0 = d0;
            this.d1 = d1;
            h = (z1-z0);
            denom = h*h;
        }
        public double eval1(double x){
            double ret = (x-z0)*(x-z1)/denom;
            ret *= ((x-z0)*d1 + (x-z1)*d0);
            return ret;
        }
        
        public double der(double x){
            double ret = (
                    (d1+d0)*(x-z0)*(x-z1) + (x-z0)*((x-z0)*d1 + (x-z1)*d0) + (x-z1)*((x-z0)*d1 + (x-z1)*d0)
                    )/denom;
            return ret;
        }
        
        public double secondDerRHS(){
            double ret = 2*(d0+2*d1)/h;
            return ret;
        }
       
        void xmin(double[] oldest, double[] upper, double[] wts, double precision){
            double x0 = oldest[0];
            double x1 = upper[0];
            double xmin = (wts[0]*x0 + wts[1]*x1)/(wts[0]+wts[1]);
            double dermin = der(xmin);
            if(Math.abs(dermin)<precision){
                oldest[0]=xmin;
                oldest[1]=dermin;
                oldest[2]=100;
                return;
            } 
            fixLimits(oldest, upper, xmin, dermin);
            xmin = (oldest[0] + upper[0])/(2);
            dermin = der(xmin);
            fixLimits(oldest, upper, xmin, dermin);
        }
        
        private void fixLimits(double[] oldest, double[] upper, double xmin, double dermin) {
            if(Math.signum(dermin) == Math.signum(oldest[1])){
                oldest[0]=xmin;
                oldest[1]=dermin;
            } else {
                upper[0]=xmin;
                upper[1]=dermin;
            }
        }
        
        abstract void estimateC(  double xmin) ;
        public abstract double eval(double x);
       
        protected double processMax( ) {
            double[] oldest = new double[]{z0, d0, 0};
            double[] upper = new double[]{z1, d1, 1};
            double[] wts = new double[]{Math.abs(d1),Math.abs(d0)};
            double xmin = z0 - d0*(z1-z0)/(d1-d0);
            double precision = 0.001*Math.abs(d0);
            for (int i = 0; i < 10; i++) {
               xmin(oldest, upper, wts, precision);
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
            fixLimits(oldest, upper, xmin, dermin);            
            wts = new double[]{Math.abs(upper[1]),Math.abs(oldest[1])};
            precision = 0.0001;
            for (int i = 0; i < 10; i++) {
               estimateC(xmin);
               xmin(oldest, upper, wts, precision);
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
            return xmin;
        }

    };

    public static class Poly4 extends Poly3{
        double C;
        double positionMax;
        double max;
        public Poly4(ZeroInfo zeroInput){
            this(zeroInput.lastZero[0],zeroInput.nextValues[0],
                    zeroInput.lastZero[1],zeroInput.nextValues[1],zeroInput.lastZero[2]);
        }
        public Poly4(double z0, double z1, double d0, double d1, double max) {
            super(z0, z1, d0, d1);
            if(d0<0){max = -max;}
            this.max= max;
            positionMax = processMax();
        }
        
        void estimateC(  double xmin) {
            double mult = (xmin-z0)*(xmin-z1);
            mult = mult*mult/denom;
            C = (max - eval1(xmin))/mult;
        } 

        public double eval(double x){
            double mult = (x-z0)*(x-z1);
            mult = mult*mult/denom;
            double ret = super.eval1(x)+C*mult;
            return ret;
        }        
        public double der(double x){
            double ret = super.der(x);
            ret += 2*C*(x-z0)*(x-z1)*(2*x-z1-z0)/denom;
            return ret;
        }
        public double secondDerRHS(){
            double ret = super.secondDerRHS();
            return ret += 2*C;        
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
        } catch (Exception e){
            throw new IllegalStateException(e);
        }
    }

    private static void readItems(   )
            throws FileNotFoundException, IOException {
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
        
        zetaCorrection1 = GSeries.correction( gSeries.basesqrtArg1);
        System.out.println( gSeries.begin + ", zetaCorrection " + zetaCorrection1);
        
        int N = Rosser.getParamInt("N");
        N = 1000200;
        int count = 0;
        zeroInput = Rosser.readZeros(baseLimit, out, zeroIn, null);
        System.out.println(Arrays.toString(zeroInput.lastZero)  +
                ", " + baseLimit + ", " + Arrays.toString(zeroInput.nextValues));
        System.arraycopy(zeroInput.nextValues, 0, lastZeroSeen1, 0, zeroIn.length);
        double[] h = new double[N];
        double[][] g = new double[N][2];
        double[] zetaMidMean = {0, 0};
        double[] zetaGramMean = {0, 0};
        int idx = 0;
        while (count < N  ) {
            int n = count + noffset;
            double upperLimit = baseLimit + (n-correction-1)* (gramIncr);
            updateZeroInput(upperLimit);
            double zeta =  getZeta(n, upperLimit, zetaGramMean);
            g[idx][0] = ((n%2==0)?-zeta:zeta);

//            System.out.println();
//            System.out.println(Arrays.toString(zeroInput.lastZero)  +
//                    ", \n" + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
//            System.out.println("gram " + zetaEst + ", " +  upperLimit + " (" + (n+1) +")");
            
            upperLimit += gramIncr/2;
            updateZeroInput(upperLimit);
            zeta = getZeta(n, upperLimit, zetaMidMean);
            h[idx] = ((n%2==0)?-zeta:zeta);
//            System.out.println("mid " + zetaEstMid + ", " +  upperLimit + " (" + (n+1) +")");
            if (count==N-1) {
                System.out.println("final n " + n );
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
        imFGramPoints(h, g);
    }

    private static double getZeta(int n, double upperLimit, double[] zetaMean) {
        double zetaEstMid = poly.eval(upperLimit);
        if(Math.abs(zeroInput.lastZero[2])>absMax){
            absMax = Math.abs(zeroInput.lastZero[2]);
            if(absMax>130){
                System.out.println();
                System.out.println(Arrays.toString(zeroInput.lastZero)  +
                        ", \n" + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
                System.out.println("secondDerRHS " + poly.secondDerRHS() + ", zetaEstMid " + zetaEstMid + " (" + (n+1) +")");
             }
        }
        double zeta = (zetaEstMid - zetaCorrection1)/2;
        zetaMean[n%2] += zeta;
        return zeta;
    }

    private static void updateZeroInput(double upperLimit) throws FileNotFoundException, IOException {
        double secondDer = Double.MIN_VALUE;
        boolean breakSeen = false;
        if(upperLimit<=zeroInput.nextValues[0]){
            zeroInput = new ZeroInfo(0, zeroInput);
        } else {
//                System.out.println("lastZeroSeen " + Arrays.toString(lastZeroSeen));
            zeroInput = Rosser.readZeros(upperLimit , out, zeroIn,  
                    zeroInput.nextValues);
            if(lastZeroSeen1[0] != zeroInput.lastZero[0]){
                // *lzs  *missedZero *current lz
                breakSeen = true;
                breaks++;
//                    System.out.println("break seen. " + Arrays.toString(lastZeroSeen1));
//                    System.out.println(Arrays.toString(zeroInput.lastZero)  +
//                            ", " + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
//                    break;
            } else if(poly != null) {
                secondDer = poly.secondDerRHS();
            }
            if(lastZeroSeen1[0]>9144.335 && lastZeroSeen1[0]<9144.374){
                /**
9144.335987630979 !=? 9144.335987630979
breakSeen  false secondDer 2.9522137235535898E17
lastZeroSeen1  [9144.335987630979, 21.916679002964564, 0.8156186116467504, 0.0]
[9144.335987630979, 21.916679002964564, 0.8156186116467504], 
9144.373942729077, [9144.449871243545, -18.66656999447713, 0.44260681775803146]
                 */
                System.out.println();
                System.out.println(lastZeroSeen1[0] + " !=? " + zeroInput.lastZero[0]);
                System.out.println("breakSeen  " + breakSeen + " secondDer " 
                        + secondDer);
                System.out.println("lastZeroSeen1  " + Arrays.toString(lastZeroSeen1));
                System.out.println(Arrays.toString(zeroInput.lastZero) + ", \n" + upperLimit + ", "
                        + Arrays.toString(zeroInput.nextValues));
                System.out.println();
            }
            System.arraycopy(zeroInput.nextValues, 0, lastZeroSeen1, 0, zeroIn.length);
            if(poly == null || secondDer == Double.MIN_VALUE){
                poly = new Poly4(zeroInput);
                return;
            }
            //poly = new Poly5(zeroInput, secondDer);
            poly = new Poly4(zeroInput);
        }
    }
    
    private static void imFGramPoints(double[] imFmid, double[][] fAtBeta ) throws IOException {
//        File imFile = new File(Rosser.getParam("conjecturesOutFile").replace(
//                "stats", "imF_Gram_"));
//        PrintStream imF_stream = new PrintStream(imFile);
//        imF_stream.println( "n, ImFGram" );
//        File reFile = new File(Rosser.getParam("conjecturesOutFile").replace(
//                "stats", "reF_Gram_"));
//        PrintStream reF_stream = new PrintStream(reFile);
//        for (int i = 0; i < fAtBeta.length; i++) {
//            reF_stream.println((i+3) + ", " + fAtBeta[i][0]);
//        }

        int correction = 0;
        if(Rosser.configParams.containsKey("correction")){
            correction = Rosser.getParamInt("correction");
        }
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        int noffset = Rosser.getParamInt("noffset");
        final BigDecimal offset = new BigDecimal(Rosser.getParam("bdoffset"));
        double begin= baseLimit + (noffset-correction-1)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        System.out.println( gSeries.begin + ", " );

        NormalizedSpline normalizedSpline = new NormalizedSpline(imFmid);
        int seriesOffset = 1, position = 1;
        normalizedSpline.evalMid(fAtBeta, seriesOffset, position);

        gSeries.rotateFtoG(fAtBeta);
        fAtBeta[0][1] =  Double.MIN_VALUE;      
        gSeries.setgAtBeta(fAtBeta);
       
        int initialPadding = 40;
        double zeta = evaluateZeta(109.99801991618585, initialPadding , gSeries);
        System.out.println( " zeta at Gram " + zeta + " cf 7.852770303334955" );
        double[] zero = {109.9434127500521, 110.10427375389713, 115.21645409737458, 
                115.35911882837084};
        double[] expectedDer = {207.28544365034014, -61.091725512779625, -7.282653909337562,
                17.960412142999786};
        for (int i = 0; i < zero.length; i++) {
            validateZero(zero[i], expectedDer[i], initialPadding, gSeries,false);
        }
        if(fAtBeta.length>=906100){
            /**
[90977.64166585186, 644.799901929005, 513.7189446414618], 
90977.65274279183, [90977.97641173516, -1487.416968799948, 23.278025237802503]
-2.6769851695451377 ** 
90977.65274279183, 14.152371466468 (905920)

 zeta at Gram 7.852770185092712 cf 7.852770303334955
 zeta -0.2391473494558046 cf 0.0
der 202.427236468793 cf 207.28544365034014
 zeta -0.01766058130254052 cf 0.0
der -70.05922387162923 cf -61.091725512779625
 zeta -7.737348062525028E-4 cf 0.0
der -7.63938703594858 cf -7.282653909337562
 zeta 0.035443910821259625 cf 0.0
der 18.39496680176594 cf 17.960412142999786

 zeta -10.065852891326445 cf 0.0
der 1198.6841000716997 cf 644.799901929005
 zeta -9.993692389363618 cf 0.0
der -1445.8444156296428 cf -1487.416968799948

             */
            double[] zero1 = {90977.64166585186, 90977.97641173516, };
            double[] expectedDer1 = {644.799901929005, -1487.416968799948, };
            for (int i = 0; i < zero1.length; i++) {
                validateZero(zero1[i], expectedDer1[i], initialPadding, gSeries,false);
            }
        }

//        for (int i = 0; i < fAtBeta.length-1; i++) {
//            imF_stream.println((i+3) + ", " + fAtBeta[i][1]);
//        }        
//        imF_stream.close();
//        reF_stream.close();
    }

    public static void main(String[] args) throws Exception{
        readItems();
    }

    public static void validateZero(double zero, double expectedDer, final int initialPadding, 
            GSeries gAtBeta, boolean checkAssert) {
        double zeta = Interpolate.evaluateZeta(zero, initialPadding, gAtBeta);
        System.out.println( " zeta " + zeta + " cf 0.0" );
        if(checkAssert){
           assertTrue(Math.abs(zeta) < 0.000001);
        }
        
        double delta = 0.01*gAtBeta.spacing;
        double zetaplus = Interpolate.evaluateZeta(zero+delta, initialPadding, gAtBeta);
        double zetaminus = Interpolate.evaluateZeta(zero-delta, initialPadding, gAtBeta);
        //244.92059950582586, 23.85164367971759, 1.0396728565623998
        double der = (zetaplus-zetaminus)/(2*delta);
        System.out.println("der " + der + " cf " + expectedDer);
        if(checkAssert){
            if(Math.abs(expectedDer-der)>0.001){
                throw new IllegalStateException("Expected " + expectedDer + " but got "
                        + der);
            }
        
        }
    }

    public static double evaluateZeta(double zero, final int initialPadding, GSeries gAtBeta) {
            double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, 
                    initialPadding, 1.6E-9, false);
    //        System.out.println("params " + zero + ", " + gAtBeta.midIdx 
    //                + Arrays.toString(gAtBeta.gAtBeta[gAtBeta.midIdx]));
            double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
            return zeta;
        }

    public static GSeries createGseries(int N) throws IOException, FileNotFoundException {
        Rosser.readConfig("data/RosserConfig.txt");
        GSeries gSeries = getRawGSeries();
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

    private static GSeries getRawGSeries() {
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
        return gSeries;
    }

   
}
