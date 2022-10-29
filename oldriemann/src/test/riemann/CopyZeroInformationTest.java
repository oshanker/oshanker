package riemann;

import java.io.*;
import java.util.Arrays;

import org.junit.Test;

import riemann.Rosser.ZeroInfo;

public class CopyZeroInformationTest {

    @Test
    public void testSkipUntil() throws FileNotFoundException, IOException {
        double[] lower = {
              244.920599505825861697, 245.916550650089922425, 246.093805026664464969,
        };
        for (int i = 0; i < lower.length; i++) {
            double[] nextValues = CopyZeroInformation.skipUntil( Interpolate.zeroIn, lower[i]);
            System.out.println(Arrays.toString(nextValues));
        }
    }

    @Test
    public void testSkipUntilBad() throws FileNotFoundException, IOException {
        double[] lower = {245.916550650089922425, 244.920599505825861697, };
        for (int i = 0; i < lower.length; i++) {
            double[] nextValues = CopyZeroInformation.skipUntil( Interpolate.zeroIn, lower[i]);
            System.out.println(Arrays.toString(nextValues));
        }
    }

    @Test
    public void testReadSingleZero() throws IOException {
        File parent = new File("out/gzeta" + Interpolate.prefix );
        if(!parent.exists()) {
        	   parent.mkdir();
        }
        BufferedReader zetaIn = null;
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
            /*
             **Gram 2 244.021159171564 1.9264980730399888
            244.021159171564-2*2*0.12179952345199391
             */
        }

        File file = new File(parent, "values.csv");
        PrintStream out = new PrintStream(file);
        ZeroInfo zeroInput = null;
        int N = 10;
        double[] nextValues =
        null
        ;

        for (int i = 0; i < N ; i++) {
            zeroInput = CopyZeroInformation.readSingleZero( Interpolate.zeroIn, nextValues);
            nextValues = zeroInput.nextValues;
            final double z0 = zeroInput.lastZero[0];
            final double z1 = zeroInput.nextValues[0];
            final double d0 = zeroInput.lastZero[1];
            final double d1 = zeroInput.nextValues[1];
            final double maxFromInput = d0>0?zeroInput.lastZero[2]:-zeroInput.lastZero[2];

            out.println("-1, " + z0 + ", 0");
            while (currentGram < z0) {
               String Zinput = zetaIn.readLine();
                System.out.println("*" + Zinput);
                String[] parsed = Zinput.split(",");
                zetaSaved = Double.parseDouble(parsed[1]);
                gramIndex = Integer.parseInt(parsed[0]);
                currentGram = baseGram + gramIndex * gramIncrement;
                System.out.println( currentGram  + ", " + zetaSaved);
            }

            Poly4 poly = new Poly4(z0,z1, d0,d1,maxFromInput);
            System.out.println(i + ", " + Arrays.toString(zeroInput.lastZero)  +
                  ", \n"  + "positionMax " + poly.positionMax 
                  + ", " + poly.eval(poly.positionMax) 
                  );
            double positionMax = poly.positionMax;
            while (currentGram < z1) {
                if (positionMax < currentGram) {
                    out.println("-100, " + positionMax + ", " + maxFromInput);
                    positionMax = Double.MAX_VALUE;
                }
                out.println( gramIndex + ", " + currentGram  + ", " + zetaSaved);
                String Zinput = zetaIn.readLine();
                System.out.println("*" + Zinput);
                String[] parsed = Zinput.split(",");
                zetaSaved = Double.parseDouble(parsed[1]);
                gramIndex = Integer.parseInt(parsed[0]);
                currentGram = baseGram + gramIndex * gramIncrement;
            }

            if(positionMax < z1) {
                out.println("-100, " + positionMax + ", " + maxFromInput);
            }
        }
        out.close();
    }

}
