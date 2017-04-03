package riemann;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class Conjectures {
	int[] signumSet;
	public int[] readItems(String filename, int sampleLength, 
			int sampleOffset, int sampleIncrement, int N)
			throws FileNotFoundException, IOException {
		/**
		 * ranked list of colos for overall cluster.
		 */
		BufferedReader zeroIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		signumSet = new int[1<<sampleLength];
//		boolean[] printed = new boolean[1<<sampleLength];
		int[] signum = new int[sampleLength];
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
				for (int i = 0; i < signum.length-1; i++) {
					signum[i] = signum[i+1];
				}
				signum[signum.length-1] = (int) ((1+Math.signum(zeta))/2);
				int checkpointForStoring = count+1-sampleLength-sampleOffset;
				if(checkpointForStoring>=0 && checkpointForStoring%sampleIncrement == 0){
					int idx = 0;
					for (int i = 0; i < signum.length; i++) {
						idx += signum[i]<<(signum.length-1-i);
					}
					signumSet[idx]++;
//					if(!printed[idx]){
//						System.out.println(idx + ":" + Arrays.toString(signum));
//						printed[idx] = true;
//					}
				}
				count++;
			}
		}
		zeroIn.close();
		return signumSet;
	}

	public static void main(String[] args) throws Exception {
		int sampleLength = 3;
		int sampleOffset = 0; 
		int sampleIncrement = 1;
		Conjectures instance = new Conjectures();
		instance.readItems("data/zeta12.csv", sampleLength, sampleOffset, sampleIncrement,
				100000);
		System.out.println(Arrays.toString(instance.signumSet));
		//[49970, 50030]
		//[39945, 10055]
		//[10025, 39975]
		//[14627, 35343, 35343, 14685]
		//[6701, 7926, 27342, 8001,  7926, 27417, 8001, 6684]
		//[1374, 5970,   730, 1951,  5909, 26692, 2012, 5361] sampleOffset = 1
		//[---,   --+,  -+-,  -++,   +--,   +-+,  ++-,   +++] 
		//[5327, 1956, 26612, 6050,  2017,   725, 5989, 1323] sampleOffset = 0


	}

}
