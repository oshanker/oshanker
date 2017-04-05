package riemann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.Arrays;

public class Conjectures {
	static String[][] descriptions = {
			{"-     ", "+     "},
			{"--   ", "-+   ", "+-   ", "++   "},
			{ "---  ", "--+  ", "-+-   ", "-++  ", "+--  ", "+-+ ", "++-  ", "+++  "}
	};
	static String[] parity = {"odd", "even"};
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
		int[] signumPoints  = new int[N];
		BufferedReader zeroIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		while (count < N) {
			String input = null;
			while ((input = zeroIn.readLine()) != null) {
				if(input == null || input.trim().length()==0){
					break;
				}
				if(input.startsWith("#") || input.startsWith("n")){
					continue;
				}
				String[] parsed = input.trim().split(",");
				double zeta = Double.parseDouble(parsed[1].trim());
				signumPoints[count] = (int) ((1+Math.signum(zeta))/2);
				count++;
			}
		}
		zeroIn.close();
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
		System.out.println(swapBits(2, 2));
		System.out.println(swapBits(5, 3));
		File file = new File("data/statsE12.csv");
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
		out.println("************** All *************");
		statistics(out, 0, 1000002);
		out.println("************** first half *************");
		statistics(out, 0, 500002);
		out.println("************** second half *************");
		statistics(out, 500000, 500002);
		out.close();
		//[ -,     +   ] 
		//[49970, 50030]
		//[39945, 10055]
		//[10025, 39975]
		//[ --,     -+,    +-,    ++ ] 
		//[14627, 35343, 35343, 14686]
		//[ 7344,  2681, 32601,  7373]
		//[ 7283, 32662,  2742,  7313]
		//[---,   --+,  -+-,  -++,   +--,   +-+,  ++-,   +++] 
		//[6701, 7926, 27342, 8001,  7926, 27417, 8001, 6684]
		//[1374, 5970,   730, 1951,  5909, 26692, 2012, 5361] even sampleOffset = 1
		//[1374, 5909,   730, 2012,  5970, 26692, 1951, 5360] reversed even
		//[5327, 1956, 26612, 6050,  2017,   725, 5989, 1323] sampleOffset = 0
		
	}

	public static void statistics(PrintStream out, int begin, int N)  {
		int sampleLength = 1;
		int sampleOffset = 0; 
		int sampleIncrement = 2;
		double[][] prob = new double[2][2];
		int maxSampleLength = 3;
		if(begin%2 != 0){
			throw new IllegalArgumentException("begin must be even");
		}
		Conjectures[][] instances = new Conjectures[maxSampleLength][2]; 
		for ( sampleLength = 1; sampleLength <= maxSampleLength; sampleLength++) {
			out.println(Arrays.toString(descriptions[sampleLength-1]));
			for (sampleOffset = 0; sampleOffset < 2; sampleOffset++) {
				Conjectures instance = new Conjectures();
				instances[sampleLength-1][sampleOffset] = instance;
				
				instance.calculateDistribution(signumPoints, 
						sampleLength, sampleOffset, sampleIncrement, N, begin);
				
				out.println(Arrays.toString(instance.signumCounts) 
						+ " " + parity[sampleOffset] + " " + instance.sampleSize);
				int length = instance.signumCounts.length;
				if(sampleOffset == 0){
					for (int i = 0; i < length; i++) {
						out.print( descriptions[sampleLength-1][length-i-1] + ", " );
					}
					out.println( " swap parity:  C1 ");
				} else {
					Conjectures instance1 = new Conjectures();
					out.println(Arrays.toString(descriptions[sampleLength-1]));
					instance1.calculateDistribution(signumPoints, 
							sampleLength, 0, 1, 1000002, 0);
					out.println(Arrays.toString(instance1.signumCounts) 
							+ " All " + instance1.sampleSize);
					for (int i = 0; i < length; i++) {
						int idx = swapBits(i, sampleLength);
						out.print( descriptions[sampleLength-1][idx] + ", " );
					}
					out.println( " All parity:  C2 ");
					
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
