package riemann;

import riemann.Rosser.ZeroInfo;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

public class CopyZeroInformation {
    /**
     * Always reads at least one line.
     * Reads one line if lower is zero or less.
     */
    public static double[] skipUntil(
          BufferedReader[] zeroIn,  double lower)
    {
        String[] input = new String[zeroIn.length];
        double[] lastValue = new double[zeroIn.length];
        double zero = Double.MAX_VALUE;
    
        try {
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
        } catch (IOException e) {
            throw new IllegalStateException("readline", e);
        }

        return lastValue;
    }
    
    /**
     * sliding window for roots
     */
    public static double[] readAndUpdateZero(
            BufferedReader[] zeroIn)
            throws IOException {
    	double[] nextValues = new double[zeroIn.length];
        String[] input = new String[zeroIn.length];
        for (int i = 0; i < input.length; i++) {
            input[i] = zeroIn[i].readLine();
        }
        for (int i = 0; i < input.length; i++) {
            try {
                if( input[i] == null) {
                	System.out.printf("done %d \n", i);
                	return null;
                }
                input[i] = input[i].trim();
                nextValues[i] = Double.parseDouble(input[i]);
            } catch (Exception e){
                System.out.println("<" + input[i] + ">");
                e.printStackTrace();
                throw e;

            }
        }  
        Rosser.update(nextValues);
		return nextValues;

    }

    public static void readAndSkip(
            BufferedReader[] zeroIn)
            throws IOException {
    	String[] input = new String[zeroIn.length];
        for (int i = 0; i < input.length; i++) {
            input[i] = zeroIn[i].readLine();
                if( input[i] == null) {
                	System.out.printf("done %d \n", i);
                	return;
                }
 
        }  

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
    	 BufferedReader[] in = Rosser.zerosFileAll("data/zerosE12.csv");
    	 
    	 int skipCount = 0; //420042;
		 for (int i = 0; i < skipCount; i++) {
    		 readAndSkip(in);
		 } 
    	 for (int i = 0; i < 3; i++) {
    		 double[] val = readAndUpdateZero(in);
    		 ZerosBuffer.put(val);
        	 showZeros();
		 } 
    }

	private static void showZeros() {
		 System.out.println(Arrays.toString(Rosser.zeros));
         System.out.println(Arrays.toString(Rosser.derivatives));
    	 System.out.println(Arrays.toString(Rosser.extrema));
    	 ZerosBuffer.printAll();
    	 System.out.println("===========");
	}

}
