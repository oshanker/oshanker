package riemann;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class SequenceWriter {

    public static void main(String[] args) throws Exception {
        String[] fileName = new String[]{"out/6KE12zeros.csv",
         "out/6KE12slopes.csv",
        "out/6KE12max.csv"};
   	    BufferedReader[] in = Rosser.zerosFileAll("data/zerosE12.csv");

        // Open the file outside the loop. 
        // Passing 'false' or omitting the second argument completely means it will OVERWRITE the file.
        try (
        		BufferedWriter writer = new BufferedWriter(new FileWriter(fileName[0], false));
        		BufferedWriter writer1 = new BufferedWriter(new FileWriter(fileName[1], false));
                BufferedWriter writer2 = new BufferedWriter(new FileWriter(fileName[2], false))
        	) {
            
       	    BufferedWriter [] writers = new BufferedWriter[] {writer,writer1,writer2};
            writers[0].write("zero_index,zero_value");
            writers[0].newLine(); 
            writers[1].write("zero_index,zero_slope");
            writers[1].newLine(); 
            writers[2].write("zero_index,zero_max");
            writers[2].newLine(); 
        	// Your sequence generation loop
            for (int seqNum = 1; seqNum <= 6000; seqNum++) {
            	double[] nextValues = CopyZeroInformation.readAndUpdateZero(in);
            	if(nextValues == null) {
            		break;
            	}
                
                // Example logic: Replace this with your actual sequence value generation
            	for (int i = 0; i < nextValues.length; i++) {
					double d = nextValues[i];
	                String line = String.format("%d, %.8f", seqNum, d);
	                writers[i].write(line);
	                writers[i].newLine(); 
				}
            }

            System.out.println("Sequence successfully written (overwritten) to " + fileName);

        } catch (IOException e) {
            System.err.println("Error writing to file: " + e.getMessage());
        }
    }
}
