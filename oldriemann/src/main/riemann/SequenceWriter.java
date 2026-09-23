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
import java.math.BigDecimal;
import java.util.Arrays;

import math.GSeries;

public class SequenceWriter {

    public static void main(String[] args) throws Exception {
        //copyZeroInfo();
    	testReadE12();
		//testReadGseries();
    }

	static void testReadGseries() throws IOException {
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
        	// calculateGSeriesE12 GSeriesTest?
	        System.out.printf("%d begin %f,  gincr %f\n", i,  begin,  gincr);
        	System.out.println(" 3rd arg " + in.readDouble());
        	for (int j = 0; j < 5; j++) {
        	System.out.println(in.readDouble());
        	}
        } else {
        	//  getSavedGSeries StaticMethods
	        System.out.printf("%d begin %f,  gincr %f 3rd arg\n", i,  begin,  gincr);
        	System.out.println(" 3rd arg " + in.readDouble());
       	    for (int j = 0; j < 5; j++) {
            	double R = in.readDouble();
    	        System.out.printf(" double R %f\n",  R);
				
			}
        }
    	System.out.println("====================");
        in.close();
        
        /*
         102566.06733986709 420045
oldriemann/src/main/riemann/Interpolate.java
Java
·
5
 (5)
import java.io.DataOutputStream;
		DataOutputStream out = outputStream( file);
    public static DataOutputStream outputStream(File file) throws FileNotFoundException {
        DataOutputStream out;
            out = new DataOutputStream(bos);


oldriemann/src/main/riemann/StaticMethods.java
Java
·
1
 (1)
        double incr = gSeries.spacing;
        double[][] gAtBeta = gSeries.gAtBeta;
         try {
             DataOutputStream out = outputStream( file);
             out.writeDouble(begin);
             out.writeDouble(incr);
             out.writeInt(gAtBeta.length);


oldriemann/src/test/math/GSeriesTest.java
Java
·
5
 (5)
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
        int k0 = 1, k1=398942;
        DataOutputStream out = null;
        File file = new File("out/" + Integer.toString(index) +"E12.dat");
Show 3 more matches


oldriemann/src/test/math/MoreGSeriesTest.java
Java
·
3
 (3)
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
        int k0 = 1, k1=398942;
        DataOutputStream out = null;
        File file = new File("out/" + Integer.toString(index) +"E12.dat");
                  out = new DataOutputStream(bos);
         */
    }

    public static DataInputStream dataInputStream(File file) throws FileNotFoundException {
        InputStream is = new FileInputStream(file);
        // create buffered input stream.
        BufferedInputStream bis = new BufferedInputStream(is);
        // create data input stream to read data in form of primitives.
        DataInputStream in = new DataInputStream(bis);
        return in;
    }

    static void testReadE12() throws Exception{
        //double t0 = gramE12[sampleIndex][0];
        int index = 102555;
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        int k0 = 1, k1=398942;
        File file = new File("data/" + Integer.toString(index) +".dat");
        InputStream is = new FileInputStream(file);
        // create buffered input stream.
        BufferedInputStream bis = new BufferedInputStream(is);

        // create data input stream to read data in form of primitives.
        DataInputStream in = new DataInputStream(bis);
        final int initialPadding = 40;
        int R = 30000+2*initialPadding;
        double begin = in.readDouble();
        double incr = in.readDouble();
        double[][] gBeta = new double[R][2];
        for (int i = 0; i < gBeta.length; i++) {
            gBeta[i][0] = in.readDouble();
            gBeta[i][1] = in.readDouble();
        }
        GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, gBeta);
        in.close();
        // line 420043
        double[] t0 = {
        		102566.06733986709, 102565.76081397697,
        		102565.95608967196,102566.06733986709
        		};
        for (double t : t0) {
            double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( 
            		t , 4, initialPadding, 1.6E-9, false);
            double zeta0 = gAtBeta.riemannZeta(gFromBLFI, t);
            System.out.println("t, " + t + " zeta0 " + zeta0);
		}

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
