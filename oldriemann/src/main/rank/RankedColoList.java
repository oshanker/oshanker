package rank;

import java.util.Arrays;
import java.util.Random;

/**
 * This class takes ranked colo lists for several sub-clusters, 
 * and generates a histogram of the Kendall Hausdorff distance to the
 * sub-cluster rankings for several possible rankings. (Generates
 * the histogram for all rankings when colo set size is 8 or less,
 * and generates the histogram for 100,000 randomly generated lists otherwise.)
 * It lists out how many of the rankings have a lower distance than
 * the overall cluster ranking 
 * @author oshanker
 *
 */
public class RankedColoList implements Comparable<RankedColoList> {

	ItemSetRankedList master;
	double kendall;
	static int count = 0;
	static int[] hist = new int[10];
	static double binWidth = 1.001/hist.length;
	/** Kendall Hausdorff distance to sub-cluster rankings, 
	 * for overall cluster ranking 
	 * */
	static double candidateDistance;
	/** How many possible rankings are better than 
	 * overall cluster ranking 
	 * */
	static int betterRanks = 0;
	
	public RankedColoList(ItemSetRankedList master, double kendall) {
		this.master = master;
		this.kendall = kendall;
	}
	
	

	public static void main(String[] args) throws Exception {
		/** read in 12 colos for and the rankings of these for the
		overall cluster, and for 7 sub-clusters */
		Item.defaultKendallHausdorfScoreForMissing = -1;
		ItemSetRankedList.initialize(12,"data/2157all.csv","data/isp.csv");
		// hardcoded, after running the List main method.
		candidateDistance = 0.1196;
		Item[] items = new Item[ItemSetRankedList.N];
		int i = 0;
		for (Item item : Item.itemMap.values()) {
			items[i++] = item;
		}
		System.out.println("items: " + Arrays.toString(items));
		if(ItemSetRankedList.N>8){
			for (int j = 0; j < 100000; j++) {
				generateRandomPerm(items);
				updateHistogram(items);
			}
		} else {
			generateAllPermutations(items.length, items);
		}
		System.out.println( "betterRanks/count: " + betterRanks + "/" + count);
		System.out.println("Histogram:" + Arrays.toString(hist));

	}
	private static <T> void swap (T[] A, int i, int j){
		T tmp = A[i];
		A[i] = A[j];
		A[j] = tmp;
	}
	public static <T> void generateRandomPerm( T[] A ) {
		Random x = new Random();
		for (int i = 0; i < A.length; i++) {
			int idx = A.length - 1 - i;
			int next = x.nextInt(idx+1);
			if(next == idx){continue;}
			swap(A, idx, next);
		}
	}
	/**
	 * generate All Permutations (Heap's algortithm)
	 * @param n
	 * @param A
	 */
	public static <T> void generateAllPermutations(int n , T[] A ) {
		if (n == 1) {
			updateHistogram(A);
			return;
		}
		for ( int i = 0; i < n - 1; i ++) {
			generateAllPermutations(n - 1, A);
			if (n%2 == 0){
				swap(A, i, n-1);
			} else {
				swap(A, 0, n-1);
			}
		}
		generateAllPermutations(n - 1, A);
	}



	private static <T> void updateHistogram(T[] A) {
		ItemSetRankedList master = new ItemSetRankedList(A.length);
		for (int i = 0; i < A.length; i++) {
			master.add((Item) A[i], i+1, -1);
		}
		Pair<Double, Double> meanKendall = ItemSetRankedList.calculateAndPrintMeanKendall(master, false);
		hist[(int) (meanKendall.first/binWidth)]++;
		if(meanKendall.first < candidateDistance){
			betterRanks++;
		}
		count++;
	}

	@Override
	public String toString() {
		return "[ " + master + ", " + kendall + "]";
	}

	public int compareTo(RankedColoList o) {
		return Double.compare(kendall, o.kendall);
	}

}
