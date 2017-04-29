package riemann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.text.NumberFormat;
import java.util.Arrays;

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
	static {
		nf.setMinimumFractionDigits(2);
		nf.setMaximumFractionDigits(2);
		nf.setGroupingUsed(false);
		intf.setMinimumIntegerDigits(6);
	}
	
	public static int swapBits(int i, int size){
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
	public static int[] reverse(int[] data) {
		int[] reversed = new int[data.length];
	    for (int i = 0; i < data.length; i++) {
			reversed[i] = data[data.length-1-i];
		}
		return reversed;
	}
	
	public static int[] readItems(String filename,  int N)
			throws FileNotFoundException, IOException {
		Rosser.hiary = true;
		int[] signumPoints  = new int[N];
		String zerosFile =Rosser.hiary?"/Users/shankero/Documents/tmp/1e12.zeros.1001_10001002":"data/zerosE12.csv";
		BufferedReader zeroIn1 = new BufferedReader(new FileReader(zerosFile));
		int count = 0;
		int signumGram = Rosser.hiary?1:0;
		double baseLimit = 244.02115917156451839965694310614387;
		//the small drifts do add up, at small values of zeta at the gram points.
		double gramIncr = 0.24359904690398668 - 1.0E-9;
		ZeroInfo zeroInput = new ZeroInfo(null,0);
		PrintStream out = null;
		while (count < N) {
			int n1 = count + (Rosser.hiary?2:1);
			double upperLimit = baseLimit + (n1-1)* (gramIncr);
			zeroInput = Rosser.readZeros(upperLimit , out, zeroIn1, zeroInput.zeroInput);
			signumPoints[count] = signumGram;
			if(zeroInput.countZeros%2 == 1){
				signumGram = signumGram==0?1:0;
			}
			count++;
		}
		zeroIn1.close();
		if(count != N){
			throw new IllegalStateException("count " + count + ", N " + N);
		}
		return signumPoints;
	}

	public int[] calculateDistribution(int[] signumPoints, int sampleLength, 
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

	public static void main(String[] args) throws Exception {
		Files.createDirectories(FileSystems.getDefault().getPath("out"));
		File file = new File("out/statsE12.csv");
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
		signumPoints = readItems("data/zetaE12.csv", 1000002);
		//PrintStream out = System.out;
		out.println("************** " + signumPoints.length + " *************");
		statistics(out, 0, 1000002);
		
	}

	public static void statistics(PrintStream out, int begin, int N)  {
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
				int parityType = Rosser.hiary?(1-sampleOffset1):sampleOffset1;
				if(sampleOffset1==1){
				out.println(Arrays.toString(descriptions[sampleLength-1]));
				}
				if(parityType == 0){
					//odd
					String counts = Arrays.toString(instance.signumCounts);
					counts = counts.replaceAll("[\\[\\]]", "");
					out.println(  parity[parityType] + " , " + counts
							+ " \\\\, "  + instance.sampleSize);
					String evencompare = " ,";
					String c1description = "";
					for (int i = 0; i < length; i++) {
						c1description += descriptions[sampleLength-1][length-i-1] + ", " ;
						int idx = swapBits(i, sampleLength);
						if(i == idx){
							evencompare += "maps~to, ";
						} else {
							evencompare += instance.signumCounts[idx] +", ";
						}
					}
					if(sampleLength == 3){
					out.println( evencompare +"\\\\ compare"+ parity[parityType] + " C2 ");
					}
					out.println(c1description + " swap parity:  C1 ");
				} else {
					Conjectures instance1 = new Conjectures();
					String counts = Arrays.toString(instance.signumCounts);
					counts = counts.replaceAll("[\\[\\]]", "");
					out.println( parity[parityType] + " , " + counts
							+ " \\\\, "  + instance.sampleSize);
					instance1.calculateDistribution(signumPoints, 
							sampleLength, 0, 1, 1000002, 0);
					String allcounts = Arrays.toString(instance1.signumCounts);
					allcounts = allcounts.replaceAll("[\\[\\]]", "");
					String evencompare = " ,";
					String compare = " ,";
					String c2description = "";
					for (int i = 0; i < length; i++) {
						int idx = swapBits(i, sampleLength);
						if(i == idx){
							evencompare += "maps~to, ";
							compare +=  "maps~to, ";
							c2description += "   self, " ;
						} else {
							evencompare += instance.signumCounts[idx] +", ";
							compare += instance1.signumCounts[idx] +", ";
							c2description += descriptions[sampleLength-1][idx] + ", " ;
						}
					}
					if(sampleLength == 3){
					out.println( evencompare +"\\\\ compare "+ parity[parityType] + "  C2 ");
					}
					out.println(" All ,   " + allcounts 
							+ "\\\\," + instance1.sampleSize);
					out.println( compare +"\\\\ compare All   C2 ");
					out.println( c2description +" All parity:  C2 ");
					
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
									*instances[1][parityIndex].sampleSize;
							out.print(nf.format( calculated) + " ");
						}
					}
					out.println(parity[parityIndex]);
				}
			}
		}
	}

}
