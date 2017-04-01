package riemann;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class Conjectures {
	public static int[] readItems(String filename, int N)
			throws FileNotFoundException, IOException {
		/**
		 * ranked list of colos for overall cluster.
		 */
		BufferedReader zeroIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		int[] signumSet = new int[8];
		int[] signum = new int[3];
		while (count < N) {
			String input = null;;
			while ((input = zeroIn.readLine()) != null) {
				if(input == null || input.trim().length()==0){
					break;
				}
				if(input.startsWith("#") || input.startsWith("n")){
					continue;
				}
				String[] parsed = input.trim().split(",");
				double zeta = Double.parseDouble(parsed[1].trim());
				signum[0] = signum[1];
				signum[1] = signum[2];
				signum[2] = (int) ((1+Math.signum(zeta))/2);
				if(count>=3){
					int idx = 4*signum[0] + 2*signum[1] + signum[2];
					signumSet[idx]++;
				}
				count++;
			}
		}
		zeroIn.close();
		return signumSet;
	}

	public static void main(String[] args) throws Exception {
		int[] signumSet = readItems("data/zeta12.csv", 100000);
		System.out.println(Arrays.toString(signumSet));
		//[6701, 7926, 27342, 8001, 7926, 27417, 8000, 6684]
	}

}
