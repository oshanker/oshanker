package rank;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;
/**
 * Represents a ranked list of items (colos).
 * The items have a string label (colo name) and a 
 * integer index. The index is used to efficiently 
 * do list operations, without the cost of
 * string comparisons.
 * The main method reads in a overall cluster ranking,
 * and several sub-cluster rankings, and calculates the
 * mean Kendall-Hausdorff distance from overall ranking
 * to the sub-cluster rankings.
 * @author oshanker
 *
 */
public class ItemSetRankedList implements RankedList{
	
	final static NumberFormat nf = NumberFormat.getInstance();
	/**
	 * number of items (colos)
	 */
	public static int N = 7;
	private static int skippedCount;
	private static double dcgMax;
	private static double dcgMin;
	int[] ranks;
	int size;
	/**
	 * add an item at rank i
	 * @param i
	 * @param rank
	 */
	public void add(Item i, int rank, double latency){
		ranks[i.index-1] = rank;
	}
	
	public void add(Iterable<String> names){
		int i = 1;
		ArrayList<Item> itemSet = new ArrayList<>(N);
		for (String name : names) {
			Item item = new Item(name);
			itemSet.add(item);
			add(item, i, -1);
			i++;
		}
	}
	
	public int rank(int index){
		return ranks[index-1];
	}
	public String toString(){
		Item[] items = new Item[size];
		for (Item item : Item.set) {
			int rank = ranks[item.index-1];
			if(rank>0){
			items[rank-1] = item;
			}
		}
		return Arrays.toString(items);
	}

	public ItemSetRankedList(int size) {
		this.size = size;
		ranks = new int[size];
	}

	/** reads in  items, and overall cluster ranking */
	public static void readItems(String filename, int N, RankedList master)
			throws FileNotFoundException, IOException {
		/**
		 * ranked list of colos for overall cluster.
		 */
		BufferedReader zeroIn = new BufferedReader(new FileReader(filename));
		int count = 0;
		while (count < N) {
			Item item = null;
			String input = null;;
			while ((input = zeroIn.readLine()) != null) {
				if(input == null || input.trim().length()==0){
					break;
				}
				if(input.startsWith("#") || input.startsWith("colo")){
					continue;
				}
				String[] parsed = input.trim().split(",");
				double latency = Integer.parseInt(parsed[6].trim());
				item = new Item(parsed[0].trim());
				master.add(item, item.index, latency);
				break;
			}
			if(item==null){break;}
			count++;
		}
		zeroIn.close();
		
	}
	/**
	 * https://en.wikipedia.org/wiki/Discounted_cumulative_gain
	 * for missing data.
	 * read in  sub-cluster cluster ranking 
	*/
	public static String readList(BufferedReader zeroIn, HashMap<String, ? super RankedList> ranksList, 
			String input, RankedList items, int N)
			throws FileNotFoundException, IOException {
		int count = 0;
		// name of sub-cluster
		String key = null;
		SUBCLUSTER:
		while (count < N) {
			Item item = null;
			if(input == null){
				input = zeroIn.readLine();
			}
			//look for a valid item belonging to the key
			double latency = -1;
			while (input  != null) {
				if( input.trim().length()==0){
					break;
				}
				if(input.startsWith("#") || input.startsWith("colo")){
					input = zeroIn.readLine();
					continue;
				}
				String[] parsed = input.trim().split(",");
				String name = parsed[0].trim();
				item = Item.itemMap.get(name);
				if(item==null){
					//System.out.println("not found " + name);
					input = zeroIn.readLine();
					continue;
				}
				latency = Integer.parseInt(parsed[1].trim());
				if(key == null){
					//initialize key only when a valid item is found
					key = parsed[3].trim();
				} else if (!key.equals(parsed[3].trim())){
					break SUBCLUSTER;
				} 
				//found key
				break;
			}
			if(item==null){break;}
			items.add(item, ++count, latency);
			input = zeroIn.readLine();
		}
		if(key != null){
			if(count==N || Item.defaultKendallHausdorfScoreForMissing > 0){
				items.setName(key);
				ranksList.put(key, items);
			} else {
//				System.out.println("bypass " + key + " found " + count);
			}
		}
		return input;
	}
	public static void main(String[] args) throws Exception {
		Item.defaultKendallHausdorfScoreForMissing = -1;
		RankedList master = initialize(12,"data/2157all.csv","data/isp.csv");
		Pair<Double, Double> mean = calculateAndPrintMeanKendall(master, true);
		System.out.println(" Mean Kendall Hausdorf dist, " +  nf.format(mean.first) 
		  + " Mean DCG, " +  nf.format(mean.second) + " skippedCount " + skippedCount + " subnets " + Item.ranksListMap.size());
	}
	static RankedList initialize(int size, String allList, String subLists) throws FileNotFoundException, IOException {
		nf.setMinimumFractionDigits(4);
		nf.setMaximumFractionDigits(4);
		nf.setGroupingUsed(false);
		N = size;
		/** reads in  items, and overall cluster ranking 
		 * Data for New York (largeaggpath 2157), 20170214
		 * Data from https://bassniumtan-hue.tan.ygrid.yahoo.com:9999/filebrowser/#/projects/amd/reports/edge/pixie_cluster_to_colo
		 * https://bassniumtan-hue.tan.ygrid.yahoo.com:9999/filebrowser/#/projects/amd/reports/edge/pixie_cluster_to_colo/20170214
		 * */
		ItemSetRankedList master = new ItemSetRankedList(N);
		readItems(allList, size, master);
		int spearmanMax = N%2==0?N*N/2:(N+1)*(N-1)/2;
		ItemSetRankedList list2 = new ItemSetRankedList(N);
		for (Item item : Item.set) {
			list2.add(item, N-item.index+1, -1);
		}
		dcgMax = DiscountedCumulativeGain.dCG(Item.set, master);	
		dcgMin = DiscountedCumulativeGain.dCG(Item.set, list2);	
		Item.kendallMax = Item.KendallHausdorf(master,list2);
		System.out.println("Master list, " + master + ", kendallMax " + Item.kendallMax);
		int subCluserCount = Integer.MAX_VALUE;
		BufferedReader zeroIn = new BufferedReader(new FileReader(subLists));
		Item.ranksListMap = new HashMap<>();
		String input = null;
		
		/** 	 read in  7 sub-cluster cluster rankings 
		 * Data for New York (largeaggpath 2157), 20170214
		 * */
		int oldCount = Item.ranksListMap.size();
		skippedCount = 0;
		for (int i = 0; i < subCluserCount; i++) {
			RankedList items = new ItemSetRankedList(N);
			input = readList(zeroIn, Item.ranksListMap, input, items, items.itemCount());
			if(oldCount == Item.ranksListMap.size() && zeroIn.ready()){
				skippedCount++;
			}
			oldCount = Item.ranksListMap.size();
			if(!zeroIn.ready()){
				System.out.println("done reading lists. skippedCount " + skippedCount);
				break;
			}
		}
		zeroIn.close();
		return master;
	}
	public static Pair<Double,Double> calculateAndPrintMeanKendall(RankedList master, boolean print) {
		File file = new File("out/out.csv");

		// if file doesnt exists, then create it
		if (!file.exists()) {
			try {
				file.createNewFile();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		PrintStream store = null;
		try {
			store = new PrintStream(file);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Set<String> subnets = Item.ranksListMap.keySet();
		int count = 0;
		double meanKendall = 0;
		double meanDCG = 0;
		int i = 0;
		if(print){
		System.out.println("subnets " + subnets.size());
		store.println(   "cluster , " + " kendall "
		+ ", " +  " DiscountedCumulativeGain");
		}
		for (String string : subnets) {
			RankedList list = (RankedList) Item.ranksListMap.get(string);
			double kendall = ((double)Item.KendallHausdorf(master,list))/Item.kendallMax;
			count++;
			meanKendall += kendall;
			double dcg = (DiscountedCumulativeGain.dCG(Item.set, list)-dcgMin)/(dcgMax-dcgMin);
			meanDCG += dcg;
			if(print){
				System.out.println("sub-cluster " + string +  ", " +  nf.format(kendall)
				  + (kendall>0.4?", **":",")+ "," + list );
				store.println( string +  ", " +  nf.format(kendall)
				+ ", " +  nf.format(dcg ));
			}
		}
		meanKendall /= count;
		meanDCG /= count;
		store.close();
		return new Pair<Double, Double>(meanKendall, meanDCG);
	}

	@Override
	public int itemCount() {
		return ranks.length;
	}

	@Override
	public double[] getValues() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getName() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void setName(String key) {
		// TODO Auto-generated method stub
		
	}

}
