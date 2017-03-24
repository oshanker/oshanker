package rank;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * Item consists of a a string label (i.e., colo name)
 * and an index for the item. The index is used to 
 * efficiently process lists of items, without the 
 * overhead of string comparisions.
 * @author oshanker
 *
 */
public class Item {
	String name;
	public int index;
	//mapping from sub-cluster name to the ranked list
	// corresponding to the sub-cluster.
	public static HashMap<String,? super RankedList> ranksListMap;
	public static double  kendallMax;
	// set of items
	public static ArrayList<Item> set = new ArrayList<>();
	// mapping from item label (colo name) to the Item object.
	static HashMap<String,Item> itemMap = new HashMap<>();
	static int cardinality;
	static double defaultKendallHausdorfScoreForMissing = -1;
	

	public Item(String name) {
		this.name = name;
		index = ++cardinality;
		itemMap.put(name, this);
		set.add(this);
	}

	/**
	 * Assume list1.rank(i) != 0 for any i
	 * @param list1: first list, full list
	 * @param list2: could be partial
	 * @return
	 */
	public static double KendallHausdorf(RankedList list1, RankedList list2){
		int mismatch = 0;
		int tie1 = 0, tie2 = 0;
		double missing = 0;
		for (int i = 1; i <= list1.itemCount(); i++) {
			for (int j = i+1; j <= list1.itemCount(); j++) {
				if(list2.rank(i)== 0 || list2.rank(j) == 0){
					//at least one element of the pair is 
					//not in list2
					if(defaultKendallHausdorfScoreForMissing > 0){
						missing += defaultKendallHausdorfScoreForMissing;
					}
					continue;
				}
				//both lists contain i and j.
				if(list1.rank(i)==list1.rank(j) && (list2.rank(i)!=list2.rank(j))){
					tie1++;
				} 
				else if(list2.rank(i)==list2.rank(j) && (list1.rank(i)!=list1.rank(j))){
					tie2++;
				} 
				else if(list1.rank(i)<list1.rank(j) && list2.rank(i)>list2.rank(j)){
					mismatch++;
				} else if(list1.rank(i)>list1.rank(j) && list2.rank(i)<list2.rank(j)){
					mismatch++;
				} 
			}
		}
		return missing + mismatch + Math.max(tie1, tie2);
		
	}

	public static int Spearman(ItemSetRankedList list1, ItemSetRankedList list2){
		int dist = 0;
		for (int i = 1; i <= list1.ranks.length; i++) {
			dist += Math.abs(list1.rank(i)-list2.rank(i));
		}
		return dist;
		
	}
	
	@Override
	public String toString() {
		return name ;
	}


	/**
	 * generate two lists, find Kendall Hausdorff distance.
	 * Used for sanity test, exercising the basic code.
	 * @param args
	 */
	public static void main(String[] args) {
		String[] names = new String[]{"A", "B", "C", "D", "E"};
		int N = names.length;
		ArrayList<Item> items = new ArrayList<>(N);
		ItemSetRankedList list1 = new ItemSetRankedList(N);
		ItemSetRankedList list2 = new ItemSetRankedList(N);
		for (int i = 0; i < N; i++) {
			Item item = new Item(names[i]);
			items.add(item);
			list1.add(item, i+1, -1);
			if(i < N-1){
				list2.add(item, N-i, -1);
			}
		}
		Item.set = items;
		System.out.println("items " + items);
		System.out.println(list1);
		System.out.println(list2);
		//int expected = N%2==0?N*N/2:(N+1)*(N-1)/2;
		System.out.println("Kendall Hausdorf dist l1 l2 " + KendallHausdorf(list1,list2));
		System.out.println("Kendall  dist l1 l1 " + KendallHausdorf(list1,list1));

	}

}
