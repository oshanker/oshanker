/**
 * 
 */
package riemann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;

/**
 * @author oshanker
 *
 */
public class Rosser {
	
	public static void readZeros(String filename,  int N, PrintStream out)
			throws FileNotFoundException, IOException {
		out.println("243.8749480149");
		BufferedReader zeroIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		String input = null;
		while (count < N && (input = zeroIn.readLine()) != null) {
			if(input == null || input.trim().length()==0){
				break;
			}
			if(input.startsWith("#") || input.startsWith("n")){
				continue;
			}
			String[] parsed = input.trim().split("\\s+");
			int n = Integer.parseInt(parsed[0].trim());
			double zeta = Double.parseDouble(parsed[1].trim());
			out.println(n +  zeta);
			count++;
		}
		zeroIn.close();
	}
	
	public static int[] readItems(String filename,  int N, PrintStream out)
			throws FileNotFoundException, IOException {
		int[] signumPoints  = new int[N];
		out.println("n-3945951431270L,good");
		BufferedReader zeroIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		String input = null;
		while (count < N && (input = zeroIn.readLine()) != null) {
			if(input == null || input.trim().length()==0){
				break;
			}
			if(input.startsWith("#") || input.startsWith("n")){
				continue;
			}
			String[] parsed = input.trim().split(",");
			int n = Integer.parseInt(parsed[0].trim());
			double zeta = Double.parseDouble(parsed[1].trim());
			if((n%2 == 1 && zeta <= 0) || (n%2 == 0 && zeta > 0)){
			   signumPoints[count] = 1;
			}
			out.println(n + "," + signumPoints[count]);
			count++;
		}
		zeroIn.close();
		return signumPoints;
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws Exception {
		File file = new File("data/zerosE12.csv");
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
		readZeros("/Users/oshanker/Google Drive/Documents/Riemann/riemann/1e12.zeros.1001_10001002.txt", 1000002, out);
		out.close();
	}

}
