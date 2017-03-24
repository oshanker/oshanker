// This class represents a gene with associated expression data
package ml.kmeans;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;

import rank.DiscountedCumulativeGain;
import rank.Item;
import rank.ItemSetRankedList;
import rank.Pair;
import rank.RankedList;

/**
 * This consists of an name (identifier), and a list of items with their 
 * values (latencies) and thier ranks.
 * 
 * @author oshanker
 *
 */
public class UserRequest implements RankedList
{
  // Data members
  private String subnetName; //gene name
  private double[] latencyValues; // Actual values of expression data
  private int[] exprRanks;     // Ranks of expression data
  
  // Constructor; sets up instance fields
  public UserRequest(String name, double[] expressionVals) 
  {
    this(name,expressionVals.length);
    for (int i = 0; i < latencyValues.length; i++){
      latencyValues[i] = expressionVals[i];
    }
    for (int i = 0; i < latencyValues.length; i++){
      for (int j = 0; j < latencyValues.length; j++){
        if (latencyValues[j] < latencyValues[i])
          exprRanks[i]++; 
      }
    }
  }

  public UserRequest(String name, int length) 
  {
    subnetName = name;
    latencyValues = new double[length];
    exprRanks = new int[length];
  }
  @Override
public String toString() {
	return " " + subnetName + " ";
}

@Override
public int hashCode() {
	final int prime = 31;
	int result = 1;
	result = prime * result + ((subnetName == null) ? 0 : subnetName.hashCode());
	return result;
}

@Override
public boolean equals(Object obj) {
	if (this == obj)
		return true;
	if (obj == null)
		return false;
	if (getClass() != obj.getClass())
		return false;
	UserRequest other = (UserRequest) obj;
	if (subnetName == null) {
		if (other.subnetName != null)
			return false;
	} else if (!subnetName.equals(other.subnetName))
		return false;
	return true;
}

// Gets expression values
  public double[] getValues() 
  {
    return latencyValues; 
  }
  
  // Gets gene name
  public String getName() 
  {
    return subnetName; 
  }
  
	@Override
	public void setName(String key) {
		subnetName = key;
	}
  
  // Computes Euclidean distance to another gene. Distance = 0 indicates
  // identical expression data, with higher values representing increasingly
  // dissimilar expression data.
  public static class euclideanDistance implements Distance
  {
	  public double distance(RankedList first, RankedList other) 
	  {
	    double dist = - 0;
	    double[] firstValues = first.getValues();
	    double[] otherValues = other.getValues();
	    for (int i = 0; i < first.itemCount(); i++) {
			double d = (firstValues[i]-otherValues[i]);
			dist += d*d;
		}
	    return Math.sqrt(dist); 
	  }
  }
  
  // Computes Spearman distance. The range of this measure
  // is between 0 and 2, with 0 representing perfectly correlated expression
  // data, 1 representing uncorrelated data, and 2 representing perfectly
  // anticorrelated data.
  public static class SpearmanDistance implements Distance{
	  public double distance(RankedList first, RankedList other) 
	  {
	    // TODO
	    // This should calculate and return the distance between the calling
	    // Gene and the argument Gene determined from the Spearman Correlation
	    return 0; 
	  }
  }

  public static class KendallDistance implements Distance{
	  public double distance(RankedList first, RankedList other) 
	  {
			double kendall = ((double)Item.KendallHausdorf(first,  other));
	        return kendall; 
	    
	  }
  }

	public int size() {
		return latencyValues.length;
	}

	@Override
	public int rank(int index) {
		return exprRanks[index-1];
	}

	@Override
	public int itemCount() {
		return exprRanks.length;
	}

	@Override
	public void add(Item i, int rank, double latency) {
		exprRanks[i.index-1] = rank;
		 latencyValues[i.index-1] = latency;
	}

	  public static void main(String[] astrArgs) throws Exception 
	  {
		  int N = 12;
		  String allList = "data/2157all.csv";
			String subLists = "data/isp.csv";
			
			RankedList master = readData(N, allList, subLists);
			Pair<Double, Double> mean = ItemSetRankedList.calculateAndPrintMeanKendall(master, true);
			System.out.println(" Mean Kendall Hausdorf dist, " +  mean.first 
			  + " Mean DCG, " +  mean.second  + " subnets " + Item.ranksListMap.size());
		    List<RankedList> members =  new ArrayList<>();
		    members.addAll((Collection<? extends RankedList>) Item.ranksListMap.values());
	  }

	public static RankedList readData(int N, String allList, String subLists)
			throws FileNotFoundException, IOException {
		RankedList master = new UserRequest("master", N);
		ItemSetRankedList.readItems(allList, N, master);
		RankedList list2 = new UserRequest("reverted", Item.set.size());
		for (Item item : Item.set) {
			list2.add(item, Item.set.size()-item.index+1, -1);
		}
		Item.kendallMax = Item.KendallHausdorf(master,list2);
		System.out.println("Master list, " + master + ", kendallMax " + Item.kendallMax);
		int subCluserCount = Integer.MAX_VALUE;
		BufferedReader zeroIn = new BufferedReader(new FileReader(subLists ));
		Item.ranksListMap = new HashMap<>();
		String input = null;

		/** 	 read in  7 sub-cluster cluster rankings 
		 * Data for New York (largeaggpath 2157), 20170214
		 * */
		int oldCount = Item.ranksListMap.size();
		int skippedCount = 0;
		for (int i = 0; i < subCluserCount; i++) {
			RankedList items = new UserRequest("list" + i, N);
			input = ItemSetRankedList.readList(zeroIn, Item.ranksListMap, input, items, items.itemCount());
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
}
