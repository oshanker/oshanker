/**
 * 
 */
package rank;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author oshanker
 *
 */
public class DiscountedCumulativeGain {
	static double defaultDiscountedCumulativeGainScoreForMissing = 0;
	static boolean DEBUG = false;
	/**
	 * relevance scores in order of ranking
	 * @param relevance
	 * @return
	 */
	static double dCG(double[] relevance){
		double dcg = 0;
		for (int i = 0; i < relevance.length; i++) {
			dcg += relevance[i]/Math.log(i+2);
		}
		dcg *= Math.log(2);
		return dcg;
	}
	
	static double dCG(Iterable<Item> set, RankedList ranks){
		int i = 0;
		double dcg = 0;
		for (Item item : set) {
			int rank = ranks.rank(item.index);
			double relevance = rank==0?defaultDiscountedCumulativeGainScoreForMissing:ranks.itemCount()+1-rank;
			dcg += relevance /Math.log(i+2);
			if(DEBUG){
				System.out.println("i " + i + " relevance " + relevance + " dcg " + dcg);
			}
			i++;
		}
		dcg *= Math.log(2);
		return dcg;
	}
	
	static double dCG(Iterable<Double> args){
		int i = 0;
		double dcg = 0;
		for (Double relevance : args) {
			dcg += relevance/Math.log(i+2);
			i++;
		}
		dcg *= Math.log(2);
		return dcg;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		List<String> masterNames = Arrays.asList("D1","D2","D3","D4","D5","D6");
		ItemSetRankedList master = new ItemSetRankedList(masterNames.size());
		master.add(masterNames);
		System.out.println(master);
		System.out.println(dCG(Item.set, master));
		ItemSetRankedList reversed = new ItemSetRankedList(masterNames.size());
		int i = reversed.size;
		for (Item item: Item.set) {
			reversed.add(item, i--, -1);
		}
		System.out.println(reversed);
		System.out.println(dCG(Item.set, reversed));
	}

}
