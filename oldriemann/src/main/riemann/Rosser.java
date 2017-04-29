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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * @author oshanker
 *
 */
public class Rosser {
	static boolean hiary = false;
	static HashMap<String, Integer> rosser = new HashMap<>();
	private static int maxS = 10;
	static class ZeroInfo{
		String zeroInput;
		int countZeros;
		public ZeroInfo(String zeroInput, int countZeros) {
			super();
			this.zeroInput = zeroInput;
			this.countZeros = countZeros;
		}
		@Override
		public String toString() {
			return "ZeroInfo [zeroInput=" + zeroInput + ", countZeros=" + countZeros + "]";
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
	
	public static ZeroInfo readZeros(  double upperLimit, PrintStream out, 
			BufferedReader zeroIn, String input)
			throws FileNotFoundException, IOException {
		ArrayList<Double> countZeros = new ArrayList<>();
		if(input == null){
			input = zeroIn.readLine();
		}
		if(input == null){
			return null;
		}
		double zero = 0;
		while (zero < upperLimit ) {
			if(input == null || input.trim().length()==0){
				System.out.println("done");
				return null;
			}
			input = input.trim();
			if(hiary){
				String[] parsed = input.split("\\s+");
				if(parsed.length<2){
					break;
				}
				zero = Integer.parseInt(parsed[0]) + Double.parseDouble(parsed[1]);
			} else {
			    zero = Double.parseDouble(input);
			}
			if(zero >= upperLimit){
				break;
			}
			print( out, zero + ",");
			countZeros.add(zero);
			input = zeroIn.readLine();
		}
		if(countZeros.size()>16){
			System.out.println("really? upperLimit " + upperLimit + " " + countZeros);
		}
		println(out, Integer.toString(countZeros.size()));
		return new ZeroInfo(input,countZeros.size());
	}
	
	public static void readItems( int N, PrintStream out)
			throws FileNotFoundException, IOException {
		String zerosFile =hiary?"/Users/oshanker/Google Drive/Documents/Riemann/riemann/1e12.zeros.1001_10001002.txt":"data/zerosE12.csv";
		//assuming that we start at a good regular odd Gram Point
		BufferedReader zeroIn = new BufferedReader(
				new FileReader(zerosFile));
		int count = 0;
		double baseLimit = 244.02115917156451839965694310614387;
		double gramIncr = 0.24359904690398668;
		ZeroInfo zeroInput = new ZeroInfo(null,0);
		boolean oldGood = true;
		boolean good = false;
		boolean inGramBlock = false;
		String interval = null;
		int signumGram = hiary?1:-1;
		int goodCount = 0, badCount = 0;
		println(out, "n-3945951431270L,good, S");
		int S = 0;//starting at regular Gram Point
		while (count < N  ) {
			int n = count + (hiary?2:1);
			double upperLimit = baseLimit + (n-1)* (gramIncr - 1.0E-9);
			if((n%2 == 1 && signumGram <= 0) || (n%2 == 0 && signumGram > 0)){
				print(out, n + ",1, " + S);
				good = true;
				goodCount++;
			} else {
				print(out, n + ",0, " + S);
				good = false;
				badCount++;
			}
			if(inGramBlock){
				interval += zeroInput.countZeros;
			}
			if(!(good ^ oldGood)){
				//both good or both bad
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
			zeroInput = readZeros(upperLimit , out, zeroIn, zeroInput.zeroInput);
			if (count==N-1) {
				System.out.println("final n " + n + " type " + good + " signumGram " + signumGram);
			}
			if (zeroInput==null) {
				break;
			}
			S += zeroInput.countZeros - 1;
			if(Math.abs(S) > maxS ){
				System.out.println("S " + S + ", n " + n + ", zeroInput.countZeros " + zeroInput.countZeros );
			}
			if(zeroInput.countZeros%2 == 1){
				signumGram = signumGram==-1?1:-1;
			}
			oldGood = good;
			count++;
		}
		double ratio = ((double)badCount)/(goodCount+badCount);
		System.out.println("goodCount " + goodCount + " badCount " + badCount + " " + ratio);
		zeroIn.close();
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws Exception {
		PrintStream out = null;
//		File file = new File("data/rosserE12.csv");
//		if (!file.exists()) {
//			try {
//				file.createNewFile();
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		}
//		try {
//			out = new PrintStream(file);
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
		
		readItems(1000002, out);
		TreeSet<String>[] stats = new TreeSet[10];
		for (int i = 0; i < stats.length; i++) {
			stats[i] = new TreeSet<String>();
		}
		
		for (String key : rosser.keySet()) {
			//String[] parsed = key.split("=");
			int idx = key.length()-2;
			if(idx < stats.length){
				stats[idx].add(key + " : " + rosser.get(key));
			} else {
				System.out.println("* " + key.length() + " " + key + " " + rosser.get(key));
			}
		}
		for (int i = 0; i < stats.length; i++) {
			System.out.println(stats[i]);
		}
//		System.out.println(rosser);
		if(out != null){out.close();}
	}

}
