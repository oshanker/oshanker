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
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import math.GSeries;
import riemann.Rosser.ZeroInfo;

public class Interpolate {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
        getZerosFile();
    }
    static BigDecimal offset;
    final static PrintStream out = null;
    static BufferedReader[] zeroIn;
    static double[] lastZeroSeen1;
    static ZeroInfo zeroInput;
    static Poly3 poly = null;
    static int breaks = 0;
    static double zetaCorrection1;
    static  double absMax = 0;
    static int correction = 0;
    static double baseLimit;
    static double gramIncr;
    static int noffset;
    private static String prefix;

    
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
            throws FileNotFoundException, IOException {
        double begin= baseLimit + (noffset-correction)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        
        zetaCorrection1 = GSeries.correction( gSeries.basesqrtArg1);
        BigDecimal tvalsi = offset.add(BigDecimal.valueOf(begin), Gram.mc);
        BigDecimal gramIndex = Gram.theta(tvalsi, Gram.mc).divide(Gram.pi, Gram.mc);
        gramIndex = 
                gramIndex.subtract(new BigDecimal("98094362213058141112271182436"), Gram.mc);
        System.out.println( gSeries.begin + ", zetaCorrection " + zetaCorrection1
                + ", gram index " + gramIndex);
        
        int N = Rosser.getParamInt("N");
        N = 41;
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

            System.out.println();
            System.out.println(Arrays.toString(zeroInput.lastZero)  +
                    ", \n" + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
            System.out.println("gram " + zeta + ", " +  upperLimit + " (" + (n+1) +")");
            
            upperLimit += gramIncr/2;
            updateZeroInput(upperLimit);
            zeta = getZeta(n, upperLimit, zetaMidMean);
            h[idx] = ((n%2==0)?-zeta:zeta);
            System.out.println("mid " + zeta + ", " +  upperLimit + " (" + (n+1) +")");
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
//        imFGramPoints(h, g);
    }

    private static double getZeta(int n, double upperLimit, double[] zetaMean) {
        double zetaEstMid = poly.eval(upperLimit);
        if(poly instanceof Poly4){
            Poly4 poly4 = (Poly4) poly;
            System.out.println("positionMax " + poly4.positionMax
                    + ", " + poly.eval(poly4.positionMax));
        }
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
        if(upperLimit<=zeroInput.nextValues[0]){
            zeroInput = new ZeroInfo(0, zeroInput);
        } else {
//                System.out.println("lastZeroSeen " + Arrays.toString(lastZeroSeen));
            zeroInput = Rosser.readZeros(upperLimit , out, zeroIn,  
                    zeroInput.nextValues);
            if(lastZeroSeen1[0] != zeroInput.lastZero[0]){
                breaks++;
//                    System.out.println("break seen. " + Arrays.toString(lastZeroSeen1));
//                    System.out.println(Arrays.toString(zeroInput.lastZero)  +
//                            ", " + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
//                    break;
            } else if(poly != null) {
            }
//            if(Math.abs(secondDer) > 100000){
//                /**
//9144.335987630979 !=? 9144.335987630979
//breakSeen  false secondDer 2.9522137235535898E17
//lastZeroSeen1  [9144.335987630979, 21.916679002964564, 0.8156186116467504, 0.0]
//[9144.335987630979, 21.916679002964564, 0.8156186116467504], 
//9144.373942729077, [9144.449871243545, -18.66656999447713, 0.44260681775803146]
//                 */
//                System.out.println();
//                System.out.println(lastZeroSeen1[0] + " !=? " + zeroInput.lastZero[0]);
//                System.out.println("breakSeen  " + breakSeen + " secondDer " 
//                        + secondDer);
//                System.out.println("lastZeroSeen1  " + Arrays.toString(lastZeroSeen1));
//                System.out.println(Arrays.toString(zeroInput.lastZero) + ", \n" + upperLimit + ", "
//                        + Arrays.toString(zeroInput.nextValues));
//                System.out.println("____________");
//            }
            System.arraycopy(zeroInput.nextValues, 0, lastZeroSeen1, 0, zeroIn.length);
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

        double begin= baseLimit + (noffset-correction-1)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        System.out.println( gSeries.begin + ", " );

        NormalizedSpline normalizedSpline = new NormalizedSpline(imFmid);
        int seriesOffset = 1, position = 1;
        normalizedSpline.evalMid(fAtBeta, seriesOffset, position);

        gSeries.rotateFtoG(fAtBeta);
        fAtBeta[0][1] =  Double.MIN_VALUE;   
        //storeG(begin, gramIncr, fAtBeta);
        //readAndValidate();

//        for (int i = 0; i < fAtBeta.length-1; i++) {
//            imF_stream.println((i+3) + ", " + fAtBeta[i][1]);
//        }        
//        imF_stream.close();
//        reF_stream.close();
    }

    public static GSeries readGSeries() throws FileNotFoundException, IOException {

        File file = new File("out/gSeries" + prefix + "/gSeries.dat");
        InputStream is = new FileInputStream(file);
        // create buffered input stream.
        BufferedInputStream bis = new BufferedInputStream(is);
        // create data input stream to read data in form of primitives.
        DataInputStream in = new DataInputStream(bis);
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
        int initialPadding = 40;
        double zeta = evaluateZeta(109.99801991618585, initialPadding , gSeries);
        System.out.println( " zeta at Gram " + zeta + " cf 7.852770303334955" );
        if(gSeries.gAtBeta.length>=906100){
/**
 
[139.51124758750822, 607.7460406187284, 345.8939561736714], 
139.54081038108032, [139.73362471736786, -7498.504547807271, 480.51147611140493]
secondDerRHS -64254.25468357052, zetaEstMid 60.46508897339184 (392)

[139.73362471736786, -7498.504547807271, 480.51147611140493], 
139.74144053703887, [139.9706475762458, 1715.3412115291642, 29.087673391667273]
secondDerRHS -146426.16864963295, zetaEstMid -58.434968749337514 (394)
 */
            double[] zero1 = {139.51124758750822, 139.73362471736786, 139.9706475762458, 90977.64166585186, 90977.97641173516, };
            double[] expectedDer1 = {607.7460406187284, -7498.504547807271, 1715.3412115291642, 644.799901929005, -1487.416968799948, };
            for (int i = 0; i < zero1.length; i++) {
                validateZero(zero1[i], expectedDer1[i], initialPadding, gSeries,false);
            }
        } else {
            double[] zero = {109.9434127500521, 110.10427375389713, 115.21645409737458, 
                    115.35911882837084};
            double[] expectedDer = {207.28544365034014, -61.091725512779625, -7.282653909337562,
                    17.960412142999786};
            /**
[109.9434127500521, 207.28544365034014, 8.283292860041835], 
110.09833499416513, [110.10427375389713, -61.091725512779625, 0.8755383742860865]
gram 0.21737417505040368, 110.09833499416513 (99)
positionMax 110.13482549953413
mid -0.3721158150659786, 110.14849253315477 (99)
positionMax 110.21670545965252

[115.21645409737458, -7.282653909337562, 0.8259045475747673], 
115.31471904908707, [115.35911882837084, 17.960412142999786, 0.38874142446558396]
gram -0.3591642240114277, 115.31471904908707 (151)
positionMax 115.39665511586186
mid 0.04981889776567724, 115.36487658807671 (151)
positionMax 115.39665511586186

[115.35911882837084, 17.960412142999786, 0.38874142446558396], 
115.41503412706635, [115.43436155012142, -17.74804957545839, 0.6247011678160491]
gram 0.14326505413359175, 115.41503412706635 (152)
positionMax 115.48193067441824
mid -0.26084094586436185, 115.46519166605599 (152)
positionMax 115.48193067441824
             */
            for (int i = 0; i < zero.length; i++) {
                validateZero(zero[i], expectedDer[i], initialPadding, gSeries,false);
            }
            
        }
    }

    private static void storeG(double begin, double incr, double[][] gAtBeta) throws IOException, FileNotFoundException {
        DataOutputStream out = null;
        // conjecturesOutFile="out/statsE28.csv"

        File file = new File("out/gSeries" + prefix + "/gSeries.dat");
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
        out.writeDouble(begin);
        out.writeDouble(incr);
        out.writeInt(gAtBeta.length);
        for (int i = 0; i < gAtBeta.length; i++) {
            out.writeDouble(gAtBeta[i][0]);
            out.writeDouble(gAtBeta[i][1]);
        }
        out.close();
    }

    public static void main(String[] args) throws Exception{
        readItems();
    }

    public static void validateZero(double zero, double expectedDer, final int initialPadding, 
            GSeries gAtBeta, boolean checkAssert) {
        double zeta = Interpolate.evaluateZeta(zero, initialPadding, gAtBeta);
        System.out.println( "# zero " + zero + " zeta " + zeta + " cf 0.0" );
        if(checkAssert){
           if(Math.abs(zeta) > 0.000001){
               throw new IllegalStateException("Expected zero but got "
                       + zeta);
           }
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
        System.out.println( gSeries.begin + ", " + zetaCorrection1);
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
        double begin= baseLimit + (noffset-correction)* (gramIncr);
        GSeries gSeries = new GSeries(1, 0, offset, begin, gramIncr);
        return gSeries;
    }

   
}
