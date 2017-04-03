package riemann;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class Conjectures {
	int[] signumCounts;
	
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
			int sampleOffset, int sampleIncrement, int N)
			throws FileNotFoundException, IOException {
		int count = 0;
		signumCounts = new int[1<<sampleLength];
//		boolean[] printed = new boolean[1<<sampleLength];
		while (count < N) {
				int checkpointForStoring = count+1-sampleLength-sampleOffset;
				if(checkpointForStoring>=0 && checkpointForStoring%sampleIncrement == 0){
					int idx = 0;
					for (int i = 0; i < sampleLength; i++) {
						idx += signumPoints[count+1-sampleLength+i]<<(sampleLength-1-i);
					}
					signumCounts[idx]++;
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
		int sampleLength = 3;
		int sampleOffset = 1; 
		int sampleIncrement = 2;
		int N = 100000;
		int[] signumPoints = readItems("data/zeta12.csv", N);
		int[] reversedPoints = reverse(signumPoints);
		Conjectures instance = new Conjectures();
		instance.calculateDistribution(signumPoints, sampleLength, sampleOffset, sampleIncrement,
				N);
		System.out.println(Arrays.toString(instance.signumCounts));
		Conjectures reversed = new Conjectures();
		reversed.calculateDistribution(reversedPoints, sampleLength, sampleOffset+1, sampleIncrement,
				N);
		System.out.println(Arrays.toString(reversed.signumCounts));
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

}
