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
import java.math.BigDecimal;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

import riemann.Rosser.GramBlock.TYPE;

/**
 * @author oshanker
 *
 */
public class Rosser {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(13);
        nf.setMaximumFractionDigits(13);
        nf.setGroupingUsed(false);
    }
	static HashMap<String, GramBlock> rosser = new HashMap<>();
    static int[][] intervalCounts = new int[2][6];
    static int[] goodIntervalCounts = new int[intervalCounts[0].length];
    static int[] badIntervalCounts = new int[intervalCounts[0].length];
	private static int maxS = 10;
	static int pEvenGood = 0, pEvenBad = 0;
	static int goodBad = 0, badGood = 0, goodGood = 0, badBad = 0;
	private static int goodCount;
	private static int badCount;
    static Map<String,String> configParams;
	
	// one instance per type of Gram Block and pattern
	public static class GramBlock{
	    enum TYPE{
	        I,II,III,NOT_REGULAR
	    }
	    TYPE type;
	    double occurrence;
	    String pattern;
        public GramBlock(String pattern, int blockZerosCount) {
            this.occurrence = 1;
            this.pattern = pattern;
            int len = pattern.length();
            if(blockZerosCount==len ){
                //we are in a regular Gram block
                if(pattern.charAt(0)=='0' && pattern.charAt(len-1)=='0'){
                    type = TYPE.III;
                }else if(pattern.charAt(0)=='0' && pattern.charAt(len-1)=='2'){
                    type = TYPE.II;
                }else if(pattern.charAt(0)=='2' && pattern.charAt(len-1)=='0'){
                    type = TYPE.I;
                }else {
                    throw new IllegalStateException();
                }
            }else{
                type = TYPE.NOT_REGULAR;
            }
        }
	    public void increment(){
	        occurrence++;
	    }
	}
	
	/**
	 * 
	 * 
	 *
	 */
	public static class ZeroInfo{
		final int countZeros;
		public final double[] lastZero;
		//zero, slope, max
		public  final double[] nextValues;
		
		private ZeroInfo(int countZeros2, double[] lastZero, double[] nextValues) {
               this.nextValues = nextValues;
               this.countZeros = countZeros2;
               this.lastZero = lastZero;
        }

        public ZeroInfo(int countZeros2, ZeroInfo zeroInfo) {
            this(countZeros2, zeroInfo.lastZero, zeroInfo.nextValues);
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
			BufferedReader[] zeroIn,  double[] nextValues)
			throws FileNotFoundException, IOException {
		ArrayList<Double> countZeros = new ArrayList<>();
		boolean skipParse = true;
		String[] input = new String[zeroIn.length];
		double[] lastValue  = new double[zeroIn.length];
		if(nextValues == null){
		    skipParse = false;
            nextValues = new double[zeroIn.length];
		    for (int i = 0; i < input.length; i++) {
	            input[i] = zeroIn[i].readLine();
            }
            if(input[0] == null || input[0].trim().length()==0){
                System.out.println("done");
                return null;
            }
		}
		double zero = 0;
		while (zero < upperLimit ) {
			if(skipParse){
			    zero = nextValues[0];
			    skipParse = false;
//	            for (int i = 0; i < input.length; i++) {
//	                lastValue[i] = nextValues[i];
//	            }
			} else {
    			input[0] = input[0].trim();
                String[] parsed = input[0].split("\\s+");
                zero = Double.parseDouble(parsed[0]);
                if(zero < 0){
                    break;
                }
    			if(parsed.length>1){
    				zero += Double.parseDouble(parsed[1]);
    			} 
                nextValues[0] = zero;
                for (int i = 1; i < input.length; i++) {
                    input[i] = input[i].trim();
                    nextValues[i] = Double.parseDouble(input[i]);
                }
			}
			if(zero >= upperLimit){
				break;
			}
			print( out, zero + ",");
            for (int i = 0; i < input.length; i++) {
                lastValue[i] = nextValues[i];
            }
			countZeros.add(zero);
            for (int i = 0; i < input.length; i++) {
                input[i] = zeroIn[i].readLine();
            }
            if(input[0] == null || input[0].trim().length()==0){
                System.out.println("done");
                return null;
            }
		}
		println(out, Integer.toString(countZeros.size()));
		return new ZeroInfo(countZeros.size(), lastValue, nextValues);
	}
	
	public static Map<String,String> readConfig(String configFile) throws IOException{
	    configParams = new HashMap<>();
        try(BufferedReader configIn = new BufferedReader(
                new FileReader(configFile))) {
            String input = configIn.readLine();
            int line = 1;
            boolean inCommentSection = false;
            while (input != null) {
                input = input.trim();
                if (input.equals("#endComment")) {
                    if (!inCommentSection) {
                        throw new IllegalStateException("not in #beginComment: line " + line);
                    }
                    inCommentSection = false;
                }
                if (input.equals("#beginComment")) {
                    if (inCommentSection) {
                        throw new IllegalStateException("nested #beginComment: line " + line);
                    }
                    inCommentSection = true;
                }
                if (input.startsWith("#") || inCommentSection || input.length() == 0) {
                    input = configIn.readLine();
                    line++;
                    continue;
                }
                String[] parsed = input.split("=");
                configParams.put(parsed[0].trim(), parsed[1].trim());
                input = configIn.readLine();
                line++;
            }
            if (inCommentSection) {
                throw new IllegalStateException("unterminated #beginComment");
            }
        }
        return configParams;
	}

    public static BufferedReader[] getZerosFile() throws FileNotFoundException {
        String zerosFile = getParam("zerosFile");
        System.out.println("zerosFile " + zerosFile);
        BufferedReader[] zeroIn = 
                {new BufferedReader(new FileReader(zerosFile))};
        return zeroIn;
    }
	
	private static void readItems(  PrintStream out, double baseLimit)
			throws FileNotFoundException, IOException {
	    BufferedReader[] zeroIn = getZerosFile();
        double gramIncr = Rosser.getParamDouble("gramIncr");
        int signumGram = Rosser.getParamInt("signumGram");
        int N = Rosser.getParamInt("N");
        int noffset = Rosser.getParamInt("noffset");
        int correction = 0;
        if(Rosser.configParams.containsKey("correction")){
            correction = Rosser.getParamInt("correction");
        }
        BigDecimal offset = null;
        if(configParams.containsKey("toffset")){
            String toffset = getParam("toffset");
            offset = new BigDecimal(toffset);
        }
        String header = getParam("header");
		//assuming that we start at a good regular (odd/even-hiary) Gram Point
		int count = 0;
		ZeroInfo zeroInput = readZeros(baseLimit, out, zeroIn, null);
		boolean oldGood = true;
		boolean good = false;
		boolean inGramBlock = false;
		String interval = null;
		goodCount = 0; badCount = 0;
		boolean evenInterval = false;
		println(out, header);
		int S = 0;//starting at regular Gram Point
		int blockZerosCount = 0;
		while (count < N  ) {
			int n = count + noffset;
			double upperLimit = baseLimit + (n-correction-1)* (gramIncr);
			if(offset != null && count%1000000 == 0){
			    double localGram = Gram.gramInterval(offset.add(BigDecimal.valueOf(upperLimit), Gram.mc));
			    System.out.println(count + ": " + upperLimit + ", " + (localGram-gramIncr));
			}
			//at entry, upperLimit=244.264758217468505 for e12
			//first zero is 244.158906912980683962
			//prev gram is 244.02115917156451839965694310614387
			//idx at entry is 3945951431270L + 2
			if((n%2 == 1 && signumGram <= 0) || (n%2 == 0 && signumGram > 0)){
				print(out, n + ", 1 , " + S + ", " + nf.format(upperLimit-gramIncr));
				good = true;
				goodCount++;
			} else {
				print(out, n + ", 0 , " + S + ", " + nf.format(upperLimit-gramIncr));
				good = false;
				badCount++;
			}
			//still dealing with old interval ( n)
			if(!(good ^ oldGood)){
				//both good or both bad
				println(out, "" );
			    if(evenInterval){
			    	throw new IllegalStateException();
			    }
				if(inGramBlock){
					interval += zeroInput.countZeros;
			        blockZerosCount += zeroInput.countZeros;
					badBad++;
				} else {
					goodGood++;
				}
			} else {
				//transition
				if(inGramBlock){
					interval += zeroInput.countZeros;
                    blockZerosCount += zeroInput.countZeros;
				    println(out,  ", exited gram Block: config " + interval);
				    if(rosser.containsKey(interval)){
				        rosser.get(interval).increment();
				    } else {
				    	rosser.put(interval, new GramBlock(interval, blockZerosCount));
				    }
				    if(oldGood || !evenInterval){
				    	throw new IllegalStateException();
				    }
				    
				    badGood++;
				} else {
					interval = Integer.toString(zeroInput.countZeros);
                    blockZerosCount = zeroInput.countZeros;
					println(out, ", entered gram Block" );
					goodBad++;
				    if(!oldGood || !evenInterval){
				    	throw new IllegalStateException();
				    }
				}
				inGramBlock = !inGramBlock;
			}
			//fetching zero count for current interval
			zeroInput = readZeros(upperLimit , out, zeroIn,  
			        zeroInput.nextValues);
			if (count==N-1) {
				System.out.println("final n " + n + " good " + good + " signumGram " + signumGram);
			}
			if (zeroInput==null) {
				break;
			}
			//store odd in 0, this is for current interval
			intervalCounts[0][zeroInput.countZeros]++;
			if(good){
			    goodIntervalCounts[zeroInput.countZeros]++;
			} else {
                badIntervalCounts[zeroInput.countZeros]++;
			}
			S += zeroInput.countZeros - 1;
			if(Math.abs(S) > maxS ){
				System.out.println("S " + S + ", n " + n + ", zeroInput.countZeros " + zeroInput.countZeros );
			}
			if(zeroInput.countZeros%2 == 1){
				signumGram = signumGram==-1?1:-1;
				evenInterval = false;
			} else {
				evenInterval = true;
			}
			oldGood = good;
			count++;
		}
		double ratio = ((double)badCount)/(goodCount+badCount);
		System.out.println("goodCount " + goodCount + " badCount " + badCount + " " + ratio);
        for (int i = 0; i < zeroIn.length; i++) {
            zeroIn[i].close();
        }
	}

    public static String getParam(String key) {
        String header = configParams.get(key);
        if(header != null) {
           header = header.substring(1, header.length()-1);
        }
        return header;
    }

    public static int getParamInt(String key) {
        int header = Integer.parseInt(configParams.get(key));
        return header;
    }

    public static double getParamDouble(String key) {
        double header = Double.parseDouble(configParams.get(key));
        return header;
    }

    private static void calculateRatio(String key,  GramBlock block, double[][] ratio) {
        int len = key.length();
        char[] chars = key.toCharArray();
        chars[0] = '2';
        chars[len-1]='0';
        String typeI = new String(chars);
        if(rosser.containsKey(typeI)){
            ratio[len-2][0]=block.occurrence;
            ratio[len-2][1] = rosser.get(new String(typeI)).occurrence;
        }
        return ;
    }

    private static void calculateRatio3(String key,  GramBlock block, double[][] ratio) {
        int len = key.length();
        int idx = key.indexOf('3');
        if(len%2==1 && idx == len/2){
            return ;
        }
        if( idx != 1 && idx != len-2){
            return ;
        }
        if(idx < len/2){
            ratio[len-2][0] += block.occurrence;
        } else {
            ratio[len-2][1] += block.occurrence;
        }
        return ;
    }
    private static String reverse(String data) {
        int length = data.length();
        char[] reversed = new char[length];
        for (int i = 0; i < length; i++) {
            reversed[length-i-1] = data.charAt(i);
        }
        return new String(reversed);
    }
    
    /**
     * @param args
     * @throws IOException 
     * @throws FileNotFoundException 
     */
    public static void main(String[] args) throws Exception {
        readConfig("data/RosserConfig.txt");
        PrintStream out = null;
        int N = Rosser.getParamInt("N");
        if (N <= 125) {
              File file = new File(getParam("conjecturesOutFile").replace("stats", "rosser"));
              if (!file.exists()) {
                  try {
                      file.createNewFile();
                  } catch (IOException e) {
                      e.printStackTrace();
                  }
              }
              try {
                  out = new PrintStream(file);
              } catch (FileNotFoundException e) {
                  e.printStackTrace();
              }
        }
        int displacementCount = 1;
        if(configParams.containsKey("displacementCount")){
            displacementCount = Integer.parseInt(configParams.get("displacementCount"));
        }
        double[][] typeIIratios = new double[10][displacementCount];
        double[][] frequencies = new double[displacementCount][intervalCounts[0].length];
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        for (int displacement = 0; displacement < displacementCount; displacement++) {
            rosser.clear();
            for (int j = 0; j < intervalCounts.length; j++) {
                for (int i = 0; i < intervalCounts[j].length; i++) {
                    intervalCounts[j][i] = 0;
                }
            }
            System.out.println("displacement " + displacement);
            readItems( out, baseLimit+(displacement-displacementCount/2)*gramIncr/10 );
            TreeSet<String>[] stats = new TreeSet[10];
            for (int i = 0; i < stats.length; i++) {
                stats[i] = new TreeSet<String>();
            }
            
            double[][]  ratio = new double[typeIIratios.length][];
            for (int i = 0; i < ratio.length; i++) {
                ratio[i] = new double[]{0, 0};
            }

            for (String key : rosser.keySet()) {
//                GramBlock reversedKey = rosser.get(reverse(key));
//                System.out.println("***** " + key + " ***** " + rosser.get(key).occurrence
//                        + ", reversed " + ((reversedKey==null)?0:reversedKey.occurrence) );
                int len = key.length();
                if(len>11){continue;}
                int idx = len -2;
                GramBlock block = rosser.get(key);
                if(block.type == TYPE.II){
                    calculateRatio(key, block, ratio);
                }
//                if(block.type == TYPE.III){
//                    calculateRatio3(key, block, ratio);
//                }
                if(idx < stats.length){
                    stats[idx].add(key + " : " + block .occurrence);
                } else {
                    System.out.println("* " + key.length() + " " + key + " " + block.occurrence);
                }
            }
            for (int i = 0; i < ratio.length; i++) {
                if(ratio[i][0]<1 || (ratio[i][1]<1) ){
                    continue;
                }
//                if( i > 5){
//                    System.out.println(stats[i] + "\n" + Arrays.toString(ratio[i]));
//                }
                typeIIratios[i][displacement] = ratio[i][0]/ratio[i][1];
            }
            /*
            for (int i = 0; i < stats.length; i++) {
                char[] typeII = new char[i+2];
                char[] typeI = new char[i+2];
                for (int j = 0; j < typeII.length; j++) {
                    if(j==0){typeII[j]='0';}
                    else if(j==typeII.length-1){typeII[j]='2';}
                    else {typeII[j]='1';}

                    if(j==0){typeI[j]='2';}
                    else if(j==typeI.length-1){typeI[j]='0';}
                    else {typeI[j]='1';}
                }
                if(rosser.containsKey(new String(typeI)) && rosser.containsKey(new String(typeII))){
//                  System.out.println(
//                          Conjectures.nf.format(((double)rosser.get(new String(typeII)))/rosser.get(new String(typeI))));
                    typeIIratios[typeII.length-2][displacement] =
                            ((double)rosser.get(new String(typeII)).occurrence)/rosser.get(new String(typeI)).occurrence;
//                    System.out.println(stats[i] );
                } else {
//                  System.out.println(stats[i] );
                }
            }
    */
     for (int i = 0; i < intervalCounts[0].length; i++) {
          System.out.print((i==0?("goodIntervalCounts &"): " &") 
                  + Conjectures.nf.format(((double)goodIntervalCounts[i])/goodCount)) ;
      }
     System.out.println(" \\\\");
      for (int i = 0; i < intervalCounts[0].length; i++) {
          System.out.print((i==0?("badIntervalCounts &"): " &") 
                  + Conjectures.nf.format(((double)badIntervalCounts[i])/badCount)) ;
      }
      System.out.println(" \\\\");
            System.out.println(" goodIntervalCounts "  + Arrays.toString(goodIntervalCounts) );
            System.out.println(" badIntervalCounts "  + Arrays.toString(badIntervalCounts) );
            for (int j = 0; j < intervalCounts.length; j++) {
                int sum = 0;
                for (int i = 0; i < intervalCounts[j].length; i++) {
                    sum += i*intervalCounts[j][i];
                }
                System.out.println(" intervalCounts " + j + " "+ Arrays.toString(intervalCounts[j]) + " sum " + sum);
                if(j==0){
                    for (int i = 0; i < intervalCounts[j].length; i++) {
                        frequencies[displacement][i] = ((double)intervalCounts[j][i])/sum;
                    }
                }
            }
            double count = goodCount+badCount;
            double pGood = ((double)goodCount)/count ;
            double pBad = ((double)badCount)/count;
            System.out.println("goodGood " + goodGood + " badGood " + badGood + " goodBad " + goodBad + " badBad " + badBad);
            System.out.println("pEvenBad " + ((double)badGood)/badCount + " pEvenGood " + ((double)goodBad)/goodCount + " ");
            System.out.println("Table 10: pEvenBad/pEvenGood " +( ((double)badGood)/badCount )/(((double)goodBad)/goodCount)
                    + " pGood/pBad " + pGood/pBad);
        }
        
        System.out.println("Table 6");

        for (int j = 0; j < typeIIratios.length; j++) {
            for (int displacement = 0; displacement < displacementCount; displacement++) {
                System.out.print((displacement==0?(j+2+" &"): "&") 
                        + Conjectures.nf.format(typeIIratios[j][displacement] ) );
            }
            System.out.println(" \\\\");
        }

        System.out.println("Product");

        int minDisplacementCount = (displacementCount+1)/2;
        for (int j = 0; j < typeIIratios.length; j++) {
            for (int displacement = minDisplacementCount; displacement < displacementCount; displacement++) {
                System.out.print((displacement==minDisplacementCount?(j+2+" &"): " &") 
                        + Conjectures.nf.format(typeIIratios[j][displacement]* typeIIratios[j][displacementCount-1-displacement]) );
            }
            System.out.println(" \\\\");
        }
//        for (int displacement = 0; displacement < displacementCount; displacement++) {
//            double fracDisplacement = ((double)(displacement-displacementCount/2))/10.0d;
//            for (int i = 0; i < intervalCounts[0].length; i++) {
//                System.out.print((i==0?(fracDisplacement+"\\delta &"): " &") 
//                        + Conjectures.nf.format(frequencies[displacement][i])) ;
//            }
//            System.out.println(" \\\\");
//        }
        //      System.out.println(rosser);
        if(out != null){out.close();}
    }

}
