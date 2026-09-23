package riemann;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;

public class SequenceWriter {

    public static void main(String[] args) throws Exception {
        //copyZeroInfo();
		String[] fileName = new String[]{"data/248714.dat",
		         "data/102555.dat",
		        "data/243.dat"};
		for (int i = 0; i < fileName.length; i++) {
			readGSeries( fileName[i], i);
		}
    }
    
    public static void readGSeries(String fileName, int i) throws IOException {
    	//
    	File file = new File(fileName);
    	System.out.println(fileName);
        DataInputStream in = dataInputStream(file);
        double begin = in.readDouble();
        double gincr = in.readDouble();
        if ( i == 0) {
	        int R = in.readInt();
	        System.out.printf("%d begin %f,  gincr %f, R %d\n", i,  begin,  gincr, R);
        	for (int j = 0; j < 5; j++) {
        	System.out.println(in.readInt());
	        System.out.printf("%d begin %f,  gincr %f\n", i,  begin,  gincr);
        	}
        } else {
        	//  getSavedGSeries StaticMethods
        	for (int j = 0; j < 5; j++) {
            	double R = in.readDouble();
    	        System.out.printf("%d begin %f,  gincr %f, double R %f\n", i,  begin,  gincr, R);
				
			}
        }
    	System.out.println("====================");
        in.close();
    }

    public static DataInputStream dataInputStream(File file) throws FileNotFoundException {
        InputStream is = new FileInputStream(file);
        // create buffered input stream.
        BufferedInputStream bis = new BufferedInputStream(is);
        // create data input stream to read data in form of primitives.
        DataInputStream in = new DataInputStream(bis);
        return in;
    }


	static void copyZeroInfo() throws FileNotFoundException {
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
