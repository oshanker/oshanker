package riemann;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import org.junit.Test;

import riemann.Interpolate.Poly4;
import riemann.Rosser.ZeroInfo;

public class CopyZeroInformationTest {

    @Test
    public void testReadSingleZero() throws FileNotFoundException, IOException {
        File file = new File("out/gzeta" + Interpolate.prefix + "/values.csv");
        PrintStream out = new PrintStream(file);
        ZeroInfo zeroInput = null;
        int N = 10;
        double[] nextValues = null;
        for (int i = 0; i < N ; i++) {
            zeroInput = copyZeroInformation.readSingleZero( Interpolate.zeroIn, nextValues);
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
