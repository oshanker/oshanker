package riemann;

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
import math.Quadratic;
import riemann.Rosser.ZeroInfo;

public class Interpolate {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }
    
    public static class Poly3{
        final double z0, z1;
        final double d0, d1;
        double denom;
        public Poly3(double z0, double z1, double d0, double d1) {
            this.z0 = z0;
            this.z1 = z1;
            this.d0 = d0;
            this.d1 = d1;
            denom = (z0-z1)*(z0-z1);
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
        
        void xmin(double[] oldest, double[] upper, double[] wts, double precision){
            double x0 = oldest[0];
            double x1 = upper[0];
            double xmin = (wts[0]*x0 + wts[1]*x1)/(wts[0]+wts[1]);
            double dermin = der(xmin);
            if(Math.abs(dermin)<precision){
                oldest[0]=xmin;
                oldest[1]=dermin;
                oldest[2]=100;
            } else if(Math.signum(dermin) == Math.signum(oldest[1])){
                oldest[0]=xmin;
                oldest[1]=dermin;
            } else {
                upper[0]=xmin;
                upper[1]=dermin;
            }
        }

        double estimateC( double max, double xmin) {
            double mult = (xmin-z0)*(xmin-z1);
            mult = mult*mult/denom;
            return (max - eval1(xmin))/mult;
        } 
    };

    public static class Poly4 extends Poly3{
        double C;
        double min;
        public Poly4(ZeroInfo zeroInput){
            this(zeroInput.lastZero[0],zeroInput.nextValues[0],zeroInput.lastZero[1],zeroInput.nextValues[1],zeroInput.lastZero[2]);
        }
        public Poly4(double z0, double z1, double d0, double d1, double max) {
            super(z0, z1, d0, d1);
            if(d0<0){max = -max;}
            double[] oldest = new double[]{z0, d0, 0};
            double[] upper = new double[]{z1, d1, 1};
            double[] wts = new double[]{Math.abs(d1),Math.abs(d0)};
            double xmin = z0 - d0*(z1-z0)/(d1-d0);
            double precision = 0.01*Math.abs(d0);
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
            C = estimateC( max, xmin);
            oldest = new double[]{z0, d0, 0};
            upper = new double[]{z1, d1, 1};
            double dermin = der(xmin);
            if(Math.signum(dermin) == Math.signum(oldest[1])){
                oldest[0]=xmin;
                oldest[1]=dermin;
            } else {
                upper[0]=xmin;
                upper[1]=dermin;
            }            
            wts = new double[]{Math.abs(upper[1]),Math.abs(oldest[1])};
            precision = 0.01;
            for (int i = 0; i < 10; i++) {
               C = estimateC( max, xmin);
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
            min = xmin;
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
    }

    public static BufferedReader[] getZerosFile() throws FileNotFoundException {
        String zerosFile = Rosser.getParam("zerosFile");
        System.out.println("zerosFile " + zerosFile);
        BufferedReader[] zeroIn = 
                {new BufferedReader(new FileReader(zerosFile)),null,null};
        String derFile = zerosFile + ".der";
        zeroIn[1] = new BufferedReader(new FileReader(derFile));
        String maxFile = zerosFile + ".max";
        zeroIn[2] = new BufferedReader(new FileReader(maxFile));
        return zeroIn;
    }

    private static void readItems(   )
            throws FileNotFoundException, IOException {
        PrintStream out = null;
        File file = new File(Rosser.getParam("conjecturesOutFile").replace(
                "stats", "imF_midGram_"));
        File reFile = new File(Rosser.getParam("conjecturesOutFile").replace(
                "stats", "reF_Gram_"));
        
        PrintStream imF_stream = new PrintStream(file);
        PrintStream reF_stream = new PrintStream(reFile);

        BufferedReader[] zeroIn = getZerosFile();
        int correction = 0;
        if(Rosser.configParams.containsKey("correction")){
            correction = Rosser.getParamInt("correction");
        }
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        int k0 = 1, k1=0;
        int noffset = Rosser.getParamInt("noffset");
        final BigDecimal offset = new BigDecimal(Rosser.getParam("bdoffset"));
        double begin= baseLimit + (noffset-correction)* (gramIncr);
        
        BigDecimal tBaseBD = BigDecimal.valueOf(begin).add(
                offset, Gram.mc);
        System.out.println("**** " + gramIncr + " cf " + Gram.gramInterval(tBaseBD));
        GSeries gSeries = new GSeries(k0, k1, offset, begin, gramIncr);
        double zetaCorrection = GSeries.correction( gSeries.basesqrtArg1);

        System.out.println( gSeries.begin + ", " + zetaCorrection);
        
        int N = Rosser.getParamInt("N");
        N = 100000;
        tBaseBD = tBaseBD.add(BigDecimal.valueOf(N*gramIncr));
        System.out.println("**** " + gramIncr + " cf " + Gram.gramInterval(tBaseBD));
        int count = 0;
        Poly4 poly = null;
        ZeroInfo zeroInput = Rosser.readZeros(baseLimit, out, zeroIn, null);
        System.out.println(Arrays.toString(zeroInput.lastZero)  +
                ", " + baseLimit + ", " + Arrays.toString(zeroInput.nextValues));
        
        double[] h = new double[N];
        double zetaMidMeanOdd = 0;
        double zetaGramMeanOdd = 0;
        double zetaMidMeanEven = 0;
        double zetaGramMeanEven = 0;
        double absMax = 0;
        int idx = 0;
        while (count < N  ) {
            int n = count + noffset;
            double upperLimit = baseLimit + (n-correction-1)* (gramIncr);
            if(upperLimit<=zeroInput.nextValues[0]){
                zeroInput = new ZeroInfo(0, zeroInput);
            } else {
                zeroInput = Rosser.readZeros(upperLimit , out, zeroIn,  
                        zeroInput.nextValues);
                poly = new Poly4(zeroInput);
            }
            if (zeroInput==null) {
                break;
            }
            double zetaEst = poly.eval(upperLimit);
            //apply correction to get F
            //what about factor of 2?
            double zeta = (zetaEst - correction)/2;
            if((n%2==0)){
                zetaGramMeanOdd += zeta;
            } else {
                zetaGramMeanEven += zeta;
            }
            reF_stream.println((n+1) + ", " + ((n%2==0)?-zeta:zeta));

//            System.out.println();
//            System.out.println(Arrays.toString(zeroInput.lastZero)  +
//                    ", " + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
//            System.out.println(zetaEst + " ** " );
            
            upperLimit += gramIncr/2;
            if(upperLimit<=zeroInput.nextValues[0]){
                zeroInput = new ZeroInfo(0, zeroInput);
            } else {
                zeroInput = Rosser.readZeros(upperLimit , out, zeroIn,  
                        zeroInput.nextValues);
                poly = new Poly4(zeroInput);
            }
            
            double zetaEstMid = poly.eval(upperLimit);
            if(Math.abs(zeroInput.lastZero[2])>absMax){
                absMax = Math.abs(zeroInput.lastZero[2]);
                if(absMax>130){
                    System.out.println();
                    System.out.println(Arrays.toString(zeroInput.lastZero)  +
                            ", " + upperLimit + ", " + Arrays.toString(zeroInput.nextValues));
                    System.out.println(zetaEst + " ** " );
                    System.out.println(upperLimit + ", " + zetaEstMid + " (" + (n+1) +")");
                 }
            }
            //apply correction to get F
            //what about factor of 2?
            zeta = (zetaEstMid - correction)/2;
            if((n%2==0)){
                zetaMidMeanOdd += zeta;
            } else {
                zetaMidMeanEven += zeta;
            }
            double imFmid = ((n%2==0)?-zeta:zeta);
            imF_stream.println((n+1) + ", " + imFmid);
            h[idx] = imFmid;
//            System.out.println(zetaEstMid + ", " +  upperLimit + " (" + (n+1) +")");
            if (count==N-1) {
                System.out.println("final n " + n );
            }
            idx++;
            count++;
        }
        System.out.println("*** zetaMidMeanOdd " + zetaMidMeanOdd/N);
        System.out.println("*** zetaGramMeanOdd " + zetaGramMeanOdd/N);
        System.out.println("*** zetaMidMeanEven " + zetaMidMeanEven/N);
        System.out.println("*** zetaGramMeanEven " + zetaGramMeanEven/N);
        double predictedSqrtArg1 = gSeries.basesqrtArg1 + gSeries.dsqrtArg1*N*gramIncr;
        zetaCorrection = GSeries.correction( predictedSqrtArg1);
        System.out.println( ": " + zetaCorrection);
        imF_stream.close();
        reF_stream.close();
        imFGramPoints(h);
    }
    
    private static void imFGramPoints(double[] imFmid) {

        NormalizedSpline normalizedSpline = new NormalizedSpline(imFmid);
        double[] actual = normalizedSpline.evalMid();
        for (int i = 0; i < 3; i++) {
            System.out.println(actual[i]);
        }
    }

    public static void main(String[] args) throws Exception{
        Rosser.readConfig("data/RosserConfig.txt");
        readItems();
    }

   
}
