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
import java.util.HashMap;

/**
 * @author oshanker
 *
 */
public class Rosser {
	static HashMap<String, Integer> rosser = new HashMap<>();
	static class ZeroInfo{
		String zeroInput;
		int countZeros;
		public ZeroInfo(String zeroInput, int countZeros) {
			super();
			this.zeroInput = zeroInput;
			this.countZeros = countZeros;
		}
		
	}
	
	static void println(PrintStream out, String message){
		if(out != null){
			out.println(message);
		}
	}
	
	static void print(PrintStream out, String message){
		if(out != null){
			out.print(message);
		}
	}
	
	public static ZeroInfo readZeros(  double upperLimit, PrintStream out, BufferedReader zeroIn, String input)
			throws FileNotFoundException, IOException {
		int countZeros = 0;
		if(input == null){
			input = zeroIn.readLine();
		}
		if(input == null){
			return new ZeroInfo(input,countZeros);
		}
		double zero = 0;
		while (zero < upperLimit ) {
			if(input == null || input.trim().length()==0){
				break;
			}
			zero = Double.parseDouble(input.trim());
			if(zero >= upperLimit){
				break;
			}
			print( out, zero + ",");
			countZeros++;
			input = zeroIn.readLine();
		}
		println(out, Integer.toString(countZeros));
		return new ZeroInfo(input,countZeros);
	}
	
	public static void readItems(String filename,  int N, PrintStream out)
			throws FileNotFoundException, IOException {
		println(out, "n-3945951431270L,good");
		BufferedReader zeroIn = new BufferedReader(new FileReader("data/zerosE12.csv"));
		BufferedReader gramIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		String input = null;
		double upperLimit = 244.02115917156451839965694310614387;
		double gramIncr = 0.24359904690398668;
		ZeroInfo zeroInput = new ZeroInfo(null,0);
		boolean oldGood = true;
		boolean good = false;
		boolean inGramBlock = false;
		String interval = null;
		while (count < N && (input = gramIn.readLine()) != null) {
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
				print(out, n + ",1" );
				good = true;
			} else {
				print(out, n + ",0" );
				good = false;
			}
			if(inGramBlock){
				interval += zeroInput.countZeros;
			}
			if(!(good ^ oldGood)){
				println(out, "" );
			} else {
				//transition
				if(inGramBlock){
				    println(out,  ", exited gram Block: config " + interval);
				    if(rosser.containsKey(interval)){
				    	rosser.put(interval, rosser.get(interval)+1);
				    } else {
				    	rosser.put(interval, 1);
				    }
				} else {
					interval = Integer.toString(zeroInput.countZeros);
					println(out, ", entered gram Block" );
				}
				inGramBlock = !inGramBlock;
			}
			zeroInput = readZeros(upperLimit, out, zeroIn, zeroInput.zeroInput);
			oldGood = good;
			upperLimit += gramIncr;
			count++;
		}
		gramIn.close();
		zeroIn.close();
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws Exception {
//		File file = new File("data/rosserE12.csv");
//		if (!file.exists()) {
//			try {
//				file.createNewFile();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}
//		PrintStream out = null;
//		try {
//			out = new PrintStream(file);
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
		readItems("data/zetaE12.csv", 1002, null);
		System.out.println(rosser);
//		out.close();
	}

}
