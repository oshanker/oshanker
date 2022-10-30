package riemann;

import java.io.*;
import java.util.Arrays;

public class createLocalData {
   public static void main(String[] args) throws Exception {
      //use 3754 as begin, 1003855 as end
      BufferedReader[] zeroIn = new BufferedReader[5];
      String[] input = new String[zeroIn.length];
      double[] nextValues = CopyZeroInformation.skipUntil(Interpolate.zeroIn, 244.15890691298068396);
      System.out.println(Arrays.toString(nextValues));
      zeroIn[0] = Interpolate.zeroIn[0];
      zeroIn[3] = Interpolate.zeroIn[1];
      zeroIn[4] = Interpolate.zeroIn[2];
      String zerosFile = "../../1e12.zeros.1001_10001002";
      String derFile = zerosFile + ".der";
      zeroIn[1] = new BufferedReader(new FileReader(derFile));
      zeroIn[1].readLine();
      String maxFile = zerosFile + ".max";
      zeroIn[2] = new BufferedReader(new FileReader(maxFile));
      zeroIn[2].readLine();
      nextValues = CopyZeroInformation.skipUntil(zeroIn, 244.367502584863394599);
      System.out.println(Arrays.toString(nextValues));
      for (int j = 0; j < 1000000; j++) {
         for (int i = 0; i < input.length; i++) {
            input[i] = zeroIn[i].readLine();
         }
      }
      //hit the bottom
      double zero = CopyZeroInformation.populateNextValues(nextValues, input);
      System.out.println("zero " + zero);
      System.out.println(Arrays.toString(nextValues));
      BufferedReader[] zeroIn3 = new BufferedReader[3];
      zeroIn3[0] = zeroIn[0];
      zeroIn3[1] = zeroIn[1];
      zeroIn3[2] = zeroIn[2];
      closeFile(zeroIn[3]);
      closeFile(zeroIn[4]);
      String baseZerosFile = Rosser.getParam("zerosFile");

      PrintStream writeToDer = new PrintStream(
            new FileOutputStream(baseZerosFile + ".der", true));
      PrintStream writeToMax = new PrintStream(
            new FileOutputStream(baseZerosFile + ".max", true));

      int additional = 0;
      for (int kk = 0; kk < additional; kk++) {
         nextValues = CopyZeroInformation.skipUntil(zeroIn3, nextValues[0]);

         if(additional-kk == 25){
            System.out.println(Arrays.toString(nextValues));
         }
//         writeToDer.println(nextValues[1]);
//         writeToMax.println(nextValues[2]);

      }

      writeToDer.close();
      writeToMax.close();
      System.out.println(Arrays.toString(nextValues));
      /*
20:46 $ wc -l zerosE12.csv
 1000003 zerosE12.csv
21:12 $ wc -l zerosE12.csv.der
    3003 zerosE12.csv.der
21:17 $ wc -l zerosE12.csv.max
    3003 zerosE12.csv.max
       */

   }

   private static void closeFile(BufferedReader zeroIn) throws IOException {
      if(zeroIn.readLine() != null) {
         throw new IllegalStateException("file not in sync");
      }
      zeroIn.close();
   }
}
