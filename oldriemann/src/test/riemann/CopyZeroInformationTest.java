package riemann;

import java.io.*;
import java.util.Arrays;

import org.junit.Test;

import riemann.Rosser.ZeroInfo;

public class CopyZeroInformationTest {

    @Test
    public void testReadSingleZero() throws FileNotFoundException, IOException {
        File parent = new File("out/gzeta" + Interpolate.prefix );
        if(!parent.exists()) {
        	   parent.mkdir();
        }
        BufferedReader zetaIn = null;
        double baseGram = 244.021159171564 - 2*2*0.12179952345199391;
        double currentGram = Double.MAX_VALUE;
        if(Interpolate.prefix.equals("E12")) {
            String zetaFile = "data/zetaE12.csv";
            zetaIn = new BufferedReader(new FileReader(zetaFile));
            String Zinput = zetaIn.readLine();
            /*
             **Gram 2 244.021159171564 1.9264980730399888
            244.021159171564-2*2*0.12179952345199391
             */
            for (int i = 0; i < 2; i++) {
                Zinput = zetaIn.readLine();
                System.out.println("*" + Zinput);
                String[] parsed = Zinput.split(",");
                double zetaSaved = Double.parseDouble(parsed[1]);
                int gramIndex = Integer.parseInt(parsed[0]);
                //baseGram
                currentGram = baseGram + gramIndex * 2 * 0.12179952345199391;
                System.out.println( currentGram  + ", " + zetaSaved);
            }
        }

        File file = new File(parent, "values.csv");
        PrintStream out = new PrintStream(file);
        ZeroInfo zeroInput = null;
        int N = 10;
        double[] nextValues = null;
        for (int i = 0; i < N ; i++) {
            zeroInput = CopyZeroInformation.readSingleZero( Interpolate.zeroIn, nextValues);
            nextValues = zeroInput.nextValues;
            final double z0 = zeroInput.lastZero[0];
            final double z1 = zeroInput.nextValues[0];
            final double d0 = zeroInput.lastZero[1];
            final double d1 = zeroInput.nextValues[1];
            final double maxFromInput = d0>0?zeroInput.lastZero[2]:-zeroInput.lastZero[2];
            Poly4 poly = new Poly4(z0,z1, d0,d1,maxFromInput);
            System.out.println(i + ", " + Arrays.toString(zeroInput.lastZero)  +
                  ", \n"  + "positionMax " + poly.positionMax 
                  + ", " + poly.eval(poly.positionMax) 
                  );
            out.println(z0 + ", 0");
            out.println(poly.positionMax + ", " + maxFromInput);
        }
        out.close();
    }

}
