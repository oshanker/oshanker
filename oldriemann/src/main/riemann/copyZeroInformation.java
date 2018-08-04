package riemann;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import riemann.Interpolate.Poly4;
import riemann.Rosser.ZeroInfo;

public class copyZeroInformation {
    
    public static ZeroInfo readZeros(   
            BufferedReader[] zeroIn,  double[] nextValues)
            throws FileNotFoundException, IOException {
        String[] input = new String[zeroIn.length];
        double[] lastValue  = new double[zeroIn.length];
        lastValue[0] = Double.NEGATIVE_INFINITY;
        double zero = 0;
        
        //populate lastValue
        if(nextValues == null){
            nextValues = new double[zeroIn.length];
            
            for (int i = 0; i < input.length; i++) {
                input[i] = zeroIn[i].readLine();
            }
            if(input[0] == null || input[0].trim().length()==0){
                System.out.println("done");
                return null;
            }
            zero = populateNextValues(lastValue, input);
            if(zero < 0 ){
                return null;
            }
        } else {
            //save last value seen
            System.arraycopy(nextValues, 0, lastValue, 0, lastValue.length);
        }
        
        //populate nextValues
        for (int i = 0; i < input.length; i++) {
            input[i] = zeroIn[i].readLine();
        }
        if(input[0] == null || input[0].trim().length()==0){
            System.out.println("done");
            return null;
        }
        zero = populateNextValues(nextValues, input);
        if(zero < 0 ){
            return null;
        }
        return new ZeroInfo(0, lastValue, nextValues);
    }


    private static double populateNextValues(double[] nextValues, String[] input) {
        double zero;
        input[0] = input[0].trim();
        String[] parsed = input[0].split("\\s+");
        zero = Double.parseDouble(parsed[0]);
        if(zero < 0){
            return zero;
        } else {
            if(parsed.length>1){
                zero += Double.parseDouble(parsed[1]);
            } 
            nextValues[0] = zero;
            for (int i = 1; i < input.length; i++) {
                input[i] = input[i].trim();
                nextValues[i] = Double.parseDouble(input[i]);
            }
        }
        return zero;
    }


    public static void main(String[] args) throws Exception {
        //use 3754 as begin
        ZeroInfo zeroInput = null;
        int N = 1010800;
        double[] nextValues = null;
        double max = 70;
        int run = 0;
        int beginRun = 0;
        for (int i = 0; i < N ; i++) {
            zeroInput = readZeros( Interpolate.zeroIn, nextValues);
            nextValues = zeroInput.nextValues;
            if(i>1003800){
                if(Math.abs(zeroInput.lastZero[1])>max){
                    run = 0;
                    beginRun = i;
                } else {
                    run++;
                }
            }
            if(run>40){
                //run at 1003813
                System.out.println("run at " + beginRun);
                break;
            }
//            if(i>=3754){
//                Poly4 poly = new Poly4(zeroInput);
//                System.out.println(i + ", " + Arrays.toString(zeroInput.lastZero)  +
//                        ", "  + "positionMax " + poly.positionMax + ", " + poly.eval(poly.positionMax) + ", "   + Arrays.toString(nextValues));
//            }
        }

    }

}
