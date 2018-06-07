package riemann;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileReader;
import org.junit.Ignore;
import org.junit.Test;

public class RosserTest {

    @Test @Ignore
    public void testGetZerosFile() {
        System.out.println("Not yet implemented");
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
