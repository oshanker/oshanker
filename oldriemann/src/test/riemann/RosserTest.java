package riemann;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

import org.junit.Ignore;
import org.junit.Test;

public class RosserTest {

    @Test @Ignore
    public void testGetZerosFile() throws IOException {
        //System.out.println("Not yet implemented");
        Rosser.readConfig("data/RosserConfig.txt");
        BufferedReader[] zeroIn = Rosser.getZerosFile();
        String input = zeroIn[0].readLine();
        double zero = Double.parseDouble(input);
        System.out.printf("first zero %f ", zero);
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        System.out.printf("gram %f \n", (baseLimit -zero)/gramIncr);
        int gramOffset = (int) Math.floor((baseLimit - zero) / gramIncr);
        double beginGram = baseLimit - gramOffset*gramIncr;
        System.out.printf("beginGram %f, check  %f\n", beginGram, beginGram + gramOffset*gramIncr);
        while (zero < beginGram) {
            input = zeroIn[0].readLine();
            zero = Double.parseDouble(input);
            System.out.printf("next zero %f \n", zero);
        }
        int count = 1; //244.158907
        double nextGram = beginGram + gramIncr;
        for (int i = 0; i < 2; i++) {
            while (zero < nextGram) {
                input = zeroIn[0].readLine();
                zero = Double.parseDouble(input);
                System.out.printf("next zero %f \n", zero);
                if(zero >= nextGram) {
                    break;
                }
                count++;
            }
            System.out.printf("i %d beginGram %f nextGram %f count %d next zero %f\n",
                    i, beginGram, nextGram, count, zero );
            // handle empty
            beginGram = nextGram;
            nextGram = beginGram + gramIncr;
            if (zero >= nextGram) {
                //while
                System.out.println("zero count!!!!");
            }
            count = 1; //more logic needed
        }


    }

    @Test
    public void testTypeIIRatios() throws Exception{
        int displacementCount = 5;
        double[][] typeIIratios = new double[10][displacementCount ];
        String ratiosFile = "data/typeIIratios.txt";
        BufferedReader ratiosIn = new BufferedReader(new FileReader(ratiosFile));
        for (int j = 0; j < typeIIratios.length; j++) {
            String val = ratiosIn.readLine();
            val = val.replaceAll("\\\\", "");
            String[] valsIn = val.split("&");
            for (int displacement = 0; displacement < displacementCount; displacement++) {
                typeIIratios[j][displacement] = Double.parseDouble(valsIn[displacement + 1]);
            }
            for (int displacement = 0; displacement < displacementCount; displacement++) {
                double val1 = Math.log(typeIIratios[j][displacement]);
                double phi = -(displacement-displacementCount/2)*Math.PI/10;
                double val2 = val1;
                if(displacement != displacementCount/2){
                    val2 = Math.tan(phi)*(1+Math.cos(2*phi)/4.5)*(j+2-0.5)/1.5;
                }
                System.out.print((displacement==0?(j+2+" & "): " & ") 
                        + Conjectures.nf.format(val1/val2 ) );
//                System.out.print((displacement==0?(j+2+" & "): " & ") 
//                        + Conjectures.nf.format(
//                                typeIIratios[j][displacement]/Math.exp(val2)) );
            }
            System.out.println(" \\\\"); 
        }
        ratiosIn.close();
    }

}
