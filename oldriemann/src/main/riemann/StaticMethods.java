package riemann;

import math.GSeries;

import java.io.*;
import java.math.BigDecimal;
import java.util.function.Function;

public class StaticMethods {
    public final static double[][] gramE12 = new double[][] {
        {243.77756012466054, 7551.727850863262},
        {7551.748966209077, 14859.592051701966},
        {14859.72037022293, 22167.406308059784},
        {22167.691772166218, 29475.511171552676},
        {29475.663172038934, 36783.50533179244},
        {36783.63456984109, 44091.59099492764},
        {44091.60596557267, 51399.45918516771},
        {51399.577359233685, 58707.39943153004},
        {58707.548750824135, 66015.36710472572},
        {66015.52014034402, 73323.37514290065},
        {73323.49152779333, 80631.16991578533},
        {80631.46291317207, 87938.90285891743},
        {87939.43429648025, 95247.28904843173},
        {95247.40567771786, 102555.26217338518},
        {102555.3770568849, 109863.22425382912},
        {109863.34843398137, 117171.2193838182},
        {117171.31980900728, 124479.11858461898},
        {124479.29118196263, 131787.26254932594},
        {131787.26255284742, 139094.915172985},
        {139095.23392166162, 146403.03320424244},
        {146403.20528840527, 153710.8534040047},
        {153711.17665307835, 161019.06253969827},
        {161019.14801568087, 168326.92904703945},
        {168327.1193762128, 175635.03465266858},
        {175635.09073467416, 182942.91104938296},
        {182943.06209106496, 190250.94773079225},
        {190251.0334453852, 197558.89402020318},
        {197559.00479763487, 204866.80626038337},
        {204866.97614781398, 212174.8229924622},
        {212174.94749592253, 219482.70295244805},
        {219482.9188419605, 226790.86141050386},
        {226790.89018592788, 234098.80158802238},
        {234098.86152782472, 241406.80651327036},
        {241406.83286765098, 248714.5783046993},
        {248714.8042054067, 256022.741415282},
    };

    static void fixLimits(double[] oldest, double[] upper, double xmin, double dermin) {
       if(Math.signum(dermin) == Math.signum(oldest[1])){
           oldest[0]=xmin;
           oldest[1]=dermin;
       } else {
           upper[0]=xmin;
           upper[1]=dermin;
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
       fixLimits(oldest, upper, xmin, dermin);
       xmin = (oldest[0] + upper[0])/(2);
       dermin = derivativeFunction.apply(xmin);
       fixLimits(oldest, upper, xmin, dermin);
   }

   public static int findFile(double t) {
      int idx = 0;
      if(t< gramE12[0][0]) {
         throw new IllegalArgumentException("out of range: " + t);
      }
      if(t> gramE12[gramE12.length-1][1]) {
         throw new IllegalArgumentException("out of range: " + t);
      }
      if(t> gramE12[gramE12.length-1][0]) {
         idx = gramE12.length-1;
      } else while(t>= gramE12[idx+1][0]) {
         idx++;
      }
      return idx;
   }

   public static GSeries getSavedGSeries(double t0, BigDecimal offset) throws FileNotFoundException, IOException {
       final int k0 = 1, k1=398942;
       final int index = (int) Math.floor(t0);
       File file = new File("data/gSeriesE12/" + Integer.toString(index) +".dat");
       InputStream is = new FileInputStream(file);
       // create buffered input stream.
       BufferedInputStream bis = new BufferedInputStream(is);
       // create data input stream to read data in form of primitives.
       DataInputStream in = new DataInputStream(bis);
       final int initialPadding = 40;
       int R = 30000+2*initialPadding;
       double begin = in.readDouble();
       double gincr = in.readDouble();
       double[][] gBeta = new double[R][2];
       for (int i = 0; i < gBeta.length; i++) {
           gBeta[i][0] = in.readDouble();
           gBeta[i][1] = in.readDouble();
       }
       GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  gincr, gBeta);
       in.close();
       return gAtBeta;
   }
}
