package riemann;

import java.io.*;
import java.util.Arrays;

import static riemann.Interpolate.zeroIn;

public class MaxInGramInterval {
   public static void getGramMax(int N) throws IOException {
      File parent = new File("out/gzeta" + Interpolate.prefix );
      if(!parent.exists()) {
         parent.mkdir();
      }
      BufferedReader zetaIn = null;
      double[] nextValues = null;
      final double gramIncrement = Interpolate.gramIncr;
      double baseGram = Interpolate.baseLimit - 2*gramIncrement;
      double currentGram = Double.MAX_VALUE;
      double zetaSaved = Double.MAX_VALUE;
      int gramIndex = Integer.MAX_VALUE;
      if(Interpolate.prefix.equals("E12")) {
         String zetaFile = "data/zetaE12.csv";
         zetaIn = new BufferedReader(new FileReader(zetaFile));
         zetaIn.readLine();
         currentGram = -1;
      }

      File file = new File(parent, "maxInGramInterval.csv");
      PrintStream out = new PrintStream(file);
      Rosser.ZeroInfo zeroInput = null;

      int gramCount = 0;
      double maxInInterval = Double.MIN_VALUE;
      for (int i = 0; i < N ; i++) {
         zeroInput = CopyZeroInformation.readSingleZero( zeroIn, nextValues);
         nextValues = zeroInput.nextValues;
         final double z0 = zeroInput.lastZero[0];
         final double z1 = zeroInput.nextValues[0];
         final double d0 = zeroInput.lastZero[1];
         final double d1 = zeroInput.nextValues[1];
         final double maxFromInput = d0>0?zeroInput.lastZero[2]:-zeroInput.lastZero[2];
         while (currentGram < z0) {
            // get next Gram
            String Zinput = zetaIn.readLine();
            String[] parsed = Zinput.split(",");
            zetaSaved = Double.parseDouble(parsed[1]);
            gramIndex = Integer.parseInt(parsed[0]);
            currentGram = baseGram + gramIndex * gramIncrement;
            if(currentGram < z0) {
               System.out.println("skipping " + currentGram + ", " + zetaSaved);
            } else {
               maxInInterval = Math.abs(zetaSaved);
            }
            if(i>0){
               throw new IllegalStateException("shouldnt be here");
            }
         }

         Poly4 poly = new Poly4(z0,z1, d0,d1,maxFromInput);
         double positionMax = poly.positionMax;
         while (currentGram < z1) {
            if (positionMax < currentGram) {
               if(maxInInterval<Math.abs(maxFromInput)){
                  maxInInterval=Math.abs(maxFromInput);
               }
               positionMax = Double.MAX_VALUE;
            }
            //write maxInInterval = zetaSaved;
            if(maxInInterval<Math.abs(zetaSaved)){
               maxInInterval=Math.abs(zetaSaved);
            }
            out.println( gramIndex + ", "  + maxInInterval);

            maxInInterval = Math.abs(zetaSaved);
            gramCount++;
            // end the gram interval
            // get next Gram
            String Zinput = zetaIn.readLine();
            String[] parsed = Zinput.split(",");
            zetaSaved = Double.parseDouble(parsed[1]);
            gramIndex = Integer.parseInt(parsed[0]);
            // incr currentGram
            currentGram = baseGram + gramIndex * gramIncrement;
         }

         if(positionMax < z1) {
            if(maxInInterval<Math.abs(maxFromInput)){
               maxInInterval=Math.abs(maxFromInput);
            }
         }
      }
      out.close();
      System.out.println("*** gramCount " + gramCount);
   }

   public static void main(String[] args) throws Exception {
      getGramMax(40);
   }
}
