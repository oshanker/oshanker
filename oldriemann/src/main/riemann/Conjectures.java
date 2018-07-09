package riemann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Map;

import riemann.Rosser.ZeroInfo;

public class Conjectures {
	static String[][] descriptions = {
			{"-     ", "+     "},
			{"--   ", "-+   ", "+-   ", "++   "},
			{ "---  ", "--+  ", "-+-   ", "-++  ", "+--  ", "+-+ ", "++-  ", "+++  "}
	};
	static String[] parity = {"Odd", "Even"};
	static int[] signumPoints;
	int[] signumCounts;
	private int sampleSize;
	static NumberFormat nf = NumberFormat.getInstance();
	static NumberFormat intf = NumberFormat.getIntegerInstance();
    private static int noffset;
	static {
		nf.setMinimumFractionDigits(6);
		nf.setMaximumFractionDigits(6);
		nf.setGroupingUsed(false);
		intf.setMinimumIntegerDigits(6);
	}
	
	private static int swapBits(int i, int size){
		int ret = 0;
		int x = i;
		int digit = 0;
		while(x>0){
			if((x&1) != 0){
				ret += 1<<(size-digit-1);
			}
			digit++;
			x >>= 1;
		}
		return ret;
	}
	private static int[] reverse(int[] data) {
		int[] reversed = new int[data.length];
	    for (int i = 0; i < data.length; i++) {
			reversed[i] = data[data.length-1-i];
		}
		return reversed;
	}
	
	private static int[] readItems(Map<String, String> configParams)
			throws FileNotFoundException, IOException {
        BufferedReader[] zeroIn1 = Rosser.getZerosFile();
        double baseLimit = Double.parseDouble(configParams.get("baseLimit"));
        double gramIncr = Double.parseDouble(configParams.get("gramIncr"));
        int signumGram = Integer.parseInt(configParams.get("signumGram"));
        signumGram = (signumGram+1)/2;
        int N = Integer.parseInt(configParams.get("N"));
        noffset = Integer.parseInt(configParams.get("noffset"));
        
		int[] signumPoints  = new int[N];
		int count = 0;
		double[] nextValues =  null;
		PrintStream out = null;
		while (count < N) {
			int n1 = count + noffset;
			double upperLimit = baseLimit + (n1-1)* (gramIncr);
			ZeroInfo zeroInput = Rosser.readZeros(upperLimit , out, zeroIn1, 
			         nextValues);
			nextValues = zeroInput.nextValues;
			signumPoints[count] = signumGram;
			if(zeroInput.countZeros%2 == 1){
				signumGram = signumGram==0?1:0;
			}
			count++;
		}
		zeroIn1[0].close();
		if(count != N){
			throw new IllegalStateException("count " + count + ", N " + N);
		}
		return signumPoints;
	}

	private int[] calculateDistribution(int[] signumPoints, int sampleLength, 
			int sampleOffset, int sampleIncrement, int N, int begin)
	{
		int count = 0;
		signumCounts = new int[1<<sampleLength];
//		boolean[] printed = new boolean[1<<sampleLength];
		while (count < N) {
			if(begin+count>=signumPoints.length){break;}
			int checkpointForStoring = count+1-sampleLength-sampleOffset;
			if(checkpointForStoring>=0 && checkpointForStoring%sampleIncrement == 0){
				int idx = 0;
				for (int i = 0; i < sampleLength; i++) {
					idx += signumPoints[begin+count+1-sampleLength+i]<<(sampleLength-1-i);
				}
				signumCounts[idx]++;
				sampleSize++;
//					if(!printed[idx]){
//						System.out.println(idx + ":" + Arrays.toString(signum));
//						printed[idx] = true;
//					}
			}
			count++;
		}
		return signumCounts;
	}

	private static void statistics(PrintStream out, int begin, int N)  {
		int sampleLength = 1;
		int sampleOffset1 = 0; 
		int sampleIncrement = 2;
		double[][] prob = new double[2][2];
		int maxSampleLength = 3;
		if(begin%2 != 0){
			throw new IllegalArgumentException("begin must be even");
		}
		Conjectures[][] instances = new Conjectures[maxSampleLength][2]; 
		for ( sampleLength = 1; sampleLength <= maxSampleLength; sampleLength++) {
			//print order of configs at top
			out.println(Arrays.toString(descriptions[sampleLength-1]));
			for (sampleOffset1 = 0; sampleOffset1 < 2; sampleOffset1++) {
				Conjectures instance = new Conjectures();
				instances[sampleLength-1][sampleOffset1] = instance;
				
				instance.calculateDistribution(signumPoints, 
						sampleLength, sampleOffset1, sampleIncrement, N, begin);
				
				int length = instance.signumCounts.length;
				// should be sampleOffset1 if !hiary
				int parityType = (1-sampleOffset1);
				if(noffset%2 == 1){
				    parityType = sampleOffset1;
				}
                String evencompare = " ,";
                double norm = instance.sampleSize;
                for (int i = 0; i < length; i++) {
//                  c1description += descriptions[sampleLength-1][length-i-1] + ", " ;
                    int idx = swapBits(i, sampleLength);
                    if(i == idx){
                        evencompare += "maps~to, ";
                    } else {
                        evencompare += nf.format(instance.signumCounts[idx]/norm) +", ";
                    }
                }
                String counts = "";
                for (int i = 0; i < instance.signumCounts.length; i++) {
                    counts += nf.format(instance.signumCounts[i]/norm) + ",";
                }
                out.println( parity[parityType] + " , " + counts
                        + " \\\\, "  + instance.sampleSize);
				if(parityType == 0){
					//odd
//					String c1description = "";
					if(sampleLength == 3){
					out.println( evencompare +"\\\\ compare "+ parity[parityType] + " C2 ");
					}
//					out.println(c1description + " swap parity:  C1 ");
				} else {
//					Conjectures instance1 = new Conjectures();
//					instance1.calculateDistribution(signumPoints, 
//							sampleLength, 0, 1, N, 0);
//					String allcounts = Arrays.toString(instance1.signumCounts);
//					allcounts = allcounts.replaceAll("[\\[\\]]", "");
					if(sampleLength == 3){
					out.println( evencompare +"\\\\ compare "+ parity[parityType] + "  C2 ");
					}
//					out.println(" All ,   " + allcounts 
//							+ "\\\\," + instance1.sampleSize);
//					out.println( compare +"\\\\ compare All   C2 ");
//					out.println( c2description +" All parity:  C2 ");
					
				}
			}
			if(sampleLength==2){
				for (int parityIndex = 0; parityIndex < parity.length; parityIndex++) {
					double sum = instances[0][parityIndex].sampleSize;
					for (int firstGram = 0; firstGram < 2; firstGram++) {
						prob[parityIndex][firstGram] = instances[0][parityIndex].signumCounts[firstGram]/sum;
					}
				}
				for (int parityIndex = 0; parityIndex < parity.length; parityIndex++) {
					
					for (int firstGram = 0; firstGram < 2; firstGram++) {
						for (int secondGram = 0; secondGram < 2; secondGram++) {
							int idx = firstGram*2 + secondGram;
							double calculated = prob[parityIndex][firstGram]* prob[(parityIndex+1)%2][secondGram]
									;
							out.print(nf.format( calculated) + " ");
						}
					}
					out.println(parity[parityIndex]);
				}
			}
		}
	}

    public static void main(String[] args) throws Exception {
        Files.createDirectories(FileSystems.getDefault().getPath("out"));
        Rosser.configParams = Rosser.readConfig("data/RosserConfig.txt");
        String conjecturesOutFile = Rosser.configParams.get("conjecturesOutFile");
        System.out.println(conjecturesOutFile);
        //get rid of quotes
        conjecturesOutFile = conjecturesOutFile.substring(1, conjecturesOutFile.length()-1);
        File file = new File(conjecturesOutFile);
        if (!file.exists()) {
            try {
                file.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        PrintStream out = null;
        try {
            out = new PrintStream(file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        signumPoints = readItems(Rosser.configParams);
        //PrintStream out = System.out;
        out.println("************** " + signumPoints.length + " *************");
        statistics(out, 0, signumPoints.length);
        
    }

}
