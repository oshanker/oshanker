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
			return null;
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
	
	public static void readItems( int N, PrintStream out)
			throws FileNotFoundException, IOException {
		//assuming that we start at a good regular odd Gram Point
		println(out, "n-3945951431270L,good");
		BufferedReader zeroIn = new BufferedReader(new FileReader("data/zerosE12.csv"));
		int count = 0;
		double baseLimit = 244.02115917156451839965694310614387;
		double gramIncr = 0.24359904690398668;
		ZeroInfo zeroInput = new ZeroInfo(null,0);
		boolean oldGood = true;
		boolean good = false;
		boolean inGramBlock = false;
		String interval = null;
		int signumGram = -1;
		while (count < N  ) {
			int n = count+1;
			if((n%2 == 1 && signumGram <= 0) || (n%2 == 0 && signumGram > 0)){
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
			double upperLimit = baseLimit + (n-1)* gramIncr;
			zeroInput = readZeros(upperLimit , out, zeroIn, zeroInput.zeroInput);
			if (zeroInput==null) {
				break;
			}
			if(zeroInput.countZeros%2 == 1){
				signumGram = signumGram==-1?1:-1;
			}
			oldGood = good;
			count++;
		}
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
		readItems(1002, null);
		System.out.println(rosser);
//		out.close();
	}

}
