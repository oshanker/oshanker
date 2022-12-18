package riemann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;

import riemann.Rosser.ZeroInfo;

public class CopyZeroInformation {
    /**
     * Always reads at least one line.
     * Reads one line if lower is zero or less.
     */
    public static double[] skipUntil(
          BufferedReader[] zeroIn,  double lower)
          throws IOException {
        String[] input = new String[zeroIn.length];
        double[] lastValue = new double[zeroIn.length];
        double zero = Double.MAX_VALUE;

        while (true) {
            for (int i = 0; i < input.length; i++) {
                input[i] = zeroIn[i].readLine();
                if(input[i] == null){
                    System.out.println("End reached");
                    return null;
                }
            }
            zero = populateNextValues(lastValue, input);
            if(zero >= lower){break;}
        }

        return lastValue;
    }

    public static ZeroInfo readSingleZero(
            BufferedReader[] zeroIn,  double[] nextValues)
            throws IOException {
        String[] input = new String[zeroIn.length];
        double[] lastValue  = new double[zeroIn.length];
        lastValue[0] = Double.NEGATIVE_INFINITY;
        double zero = 0;
        
        //populate lastValue
        if(nextValues == null){
            nextValues = new double[zeroIn.length];

            System.out.println("=====");
            for (int i = 0; i < input.length; i++) {
                input[i] = zeroIn[i].readLine();
                System.out.println(i + " " + input[i]);
            }
            System.out.println("=====");
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


    /**
     *
     * @param nextValues mutable, gets populated
     * @param input non-mutable
     * @return zero which has been read
     */
    static double populateNextValues(double[] nextValues, String[] input) {
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
                try {
                    input[i] = input[i].trim();
                    nextValues[i] = Double.parseDouble(input[i]);
                } catch (Exception e){
                    System.out.println("<" + input[i] + ">");
                    char[] ch = input[i].toCharArray();
                    for(int j = 0; j < ch.length; j++){
                        System.out.printf("char at %d index is: %d\n" , j, (int)ch[j]);
                    }
                    throw e;

                }
            }
        }
        return zero;
    }


    public static void main(String[] args) throws Exception {
        //use 3754 as begin, 1003855 as end
        File file = new File("out/gSeries" + Interpolate.prefix + "/zeros.dat");
//        DataOutputStream out = Interpolate.outputStream( file);
        ZeroInfo zeroInput = null;
        int N = 1003814;
        double[] nextValues = null;
        double max = 70;
        int run = 0;
        int beginRun = 0;
        for (int i = 0; i < N ; i++) {
            zeroInput = readSingleZero( Interpolate.zeroIn, nextValues);
            nextValues = zeroInput.nextValues;
//            if(i>1003800){
//                if(Math.abs(zeroInput.lastZero[1])>max){
//                    run = 0;
//                    beginRun = i;
//                } else {
//                    run++;
//                }
//            }
//            if(run>40){
//                //run at 1003813
//                System.out.println("run at " + beginRun);
//                break;
//            }
            if(i<3792 ){
                continue;
            }
            final double z0 = zeroInput.lastZero[0];
            final double z1 = zeroInput.nextValues[0];
            final double d0 = zeroInput.lastZero[1];
            final double d1 = zeroInput.nextValues[1];
            final double maxFromInput = d0>0?zeroInput.lastZero[2]:-zeroInput.lastZero[2];
            if(i==3792 || i==1003813){
                //gSeries.begin 476.85026008636953
                Poly4 poly = new Poly4(z0,z1, d0,d1,maxFromInput);
                System.out.println(i + ", " + Arrays.toString(zeroInput.lastZero)  +
                      ", \n"  + "positionMax " + poly.positionMax 
                      + ", " + poly.eval(poly.positionMax) 
                      + ", \n"   + Arrays.toString(nextValues));
            }
//            for (int j = 0; j < zeroInput.lastZero.length; j++) {
//                out.writeDouble(zeroInput.lastZero[j]);
//            }
//            out.writeDouble(poly.positionMax);
        }
//        out.close();
//        DataInputStream in = Interpolate.dataInputStream( file);
//        double[] tmin = new double[4];
//        for (int i = 0; i < N-3792 ; i++) 
//        {
//            for (int i1 = 0; i1 < tmin.length; i1++) 
//            {
//                tmin[i1] = in.readDouble();
//            }
//            System.out.println(Arrays.toString(tmin));
//        }
//        in.close();
    }

}
