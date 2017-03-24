// Class for computing and representing k-means clustering of expression data.
package ml.kmeans;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import rank.DiscountedCumulativeGain;
import rank.Item;
import rank.ItemSetRankedList;
import rank.Pair;
import rank.RankedList;

public class KMeans 
{
   private Distance dist;
  
  // loads expression data from file <fileName>. The genes array is
  // filled with the genes from the dataset. The clusters array is not filled;
  // to fill it, call performClustering.
  public static List<RankedList> readPSCFile(String fileName) 
  {
    BufferedReader reader;
    int numGenes;
    String[] splitLine;
    //double[] exprValues;
    List<RankedList> members =  new ArrayList<>();
    // Creates a new KMeans object by reading in all of the genes
    // and expression levels located in filename
    try 
    {
      reader = new BufferedReader(new FileReader(fileName));
      
      // Count the number of lines to determine how many genes are present
      for (numGenes = 0; reader.readLine() != null; numGenes++);
      
      // Close and then re-open the file now that we know its length
      reader.close();
      reader = new BufferedReader(new FileReader(fileName));
      
      double[] exprValues = null;
      // Now, read in each line and create the corresponding Gene object
      for (int i = 0; i < numGenes; i++) 
      {
        String line = reader.readLine();
        // The files are tab-delimited, so split on tabs (\t)
        splitLine = line.split("\t");
        
        if (exprValues == null) {
            //ret.geneExpressionSize = ;
            exprValues = new double[splitLine.length - 1];
            System.out.println("exprValues "  + exprValues.length);
        }
        // The first entry in the parts array is the gene name, the rest
        // are expression levels for that gene
        if ( (exprValues.length != (splitLine.length - 1))){
          System.out.println(exprValues.length + " != " + (splitLine.length - 1));
        }
        for (int j = 0; j < exprValues.length; j++)
          exprValues[j] = Double.parseDouble(splitLine[j + 1]);
        
        // Finally, create the Gene in the array
         members.add( new UserRequest(splitLine[0], exprValues)); 
      }
      
      // Lastly, close the file
      reader.close(); 
    }
    catch (FileNotFoundException e) 
    {
      System.out.println("ERROR:  Unable to locate " + fileName + ".");
      System.exit(0); 
    }
    catch (IOException e) 
    {
      System.out.println("ERROR:  Unable to read from " + fileName + ".");
      System.exit(0);
    }
	return members; 
  }
  
  // Perform k-means clustering with the specified number of clusters and
  // distance metric. The "metric" parameter can take two values: "euclid" for
  // Euclidean distance, or "spearman" for Spearman distance.
  public List<Cluster> performClustering(List<Cluster> clusters, List<RankedList> members ) 
  {
    double[] oldRadius = null;
    for (int i = 0; i < 20; i++) {
        //assign
		for (Cluster cluster: clusters) {
			cluster.generateCentroid();
		}
    	double[] newRadius = assign( clusters, members);
//    	 for (Cluster cluster : clusters) {
//    		System.out.println(cluster);
//    	}
    	
        //if not first, check if changed
    	 if(oldRadius != null){
    		 double change = 0;
    		 for (int j = 0; j < newRadius.length; j++) {
				change += Math.abs(newRadius[j]-oldRadius[j]);
			}
    		 if(change < 0.01){break;}
    	 }
		 oldRadius = newRadius;
	}
    
	return clusters;
  }

	public List<Cluster> initializeClusters(int numClusters, List<RankedList> members) {
		int blocksize = members.size()/(numClusters);
	    List<Cluster> clusters = new ArrayList<>();
	    //init
	    for (int i = 0; i < numClusters; i++) {
			Cluster cluster = new Cluster("C_" + numClusters + "_" + i);
			cluster.add(members.get(i*blocksize + blocksize/2));
			clusters.add(cluster);
		}
		return clusters;
	}
  
	private double tryToBreakUpClusters(List<Cluster> clusters2, double maxRadius, List<RankedList> members ){
		List<Cluster> toRemove = new ArrayList<>();
		List<Cluster> toAdd = new ArrayList<>();
		for (Cluster cluster: clusters2) {
			if(cluster.getRadius()<maxRadius){
				continue;
			}
			System.out.println("fixing " + cluster);

			double originalRadius = cluster.getRadius();
			List<Cluster> brokenUp = breakUpSpecifiedCluster(cluster);
			//TODO: need to fix re-assessment of clustering
			performClustering(  brokenUp, members);
			if(brokenUp.get(0).getRadius()< originalRadius && brokenUp.get(1).getRadius()< originalRadius ){
				//good
				System.out.println("fixed");
				toRemove.add(cluster);
				toAdd.addAll(brokenUp);
				
			} else {
				//bad
				System.out.println("not fixed");
			}
		}
		if(toRemove.size()>0){
			clusters2.removeAll(toRemove);
			clusters2.addAll(toAdd);
			performClustering(  clusters2, members);
		}
		double maxSeen = 0;
		for (Cluster cluster: clusters2) {
			if(maxSeen<cluster.getRadius()){
				maxSeen = cluster.getRadius();
			}
		}
		return maxSeen;
	}

	private static List<Cluster> breakUpSpecifiedCluster(Cluster cluster) {
		List<Cluster> brokenUp = new ArrayList<>();
		for (int i = 0; i < 2; i++) {
			Cluster cluster1 = new Cluster(cluster.getName() + "_" + i);
			RankedList member = i==0?cluster.getCentroid():cluster.getFarthest();
			cluster1.add(member);
			brokenUp.add(cluster1);
		}
		return brokenUp;
	}

  /**
   * calculates centroids, empties clusters, and assigns from set of all.
   * clusters is already populated.
   * @param metric
   * @param clusters2
   * @return
   */
	private  double[] assign( List<Cluster> clusters2,
			List<RankedList> members2) {
	    double[] newRadius = new double[clusters2.size()];
		for (Cluster cluster: clusters2) {
			cluster.clear();
		}
		for (RankedList member: members2) {
	    	Cluster min = null;
	    	double minDist = Double.MAX_VALUE;
	        for (Cluster cluster: clusters2) {
	    		double distance = dist.distance(member, cluster.getCentroid());
	    		if(distance < minDist){
	    			minDist = distance;
	    			min = cluster;
	    		}
	    	}
	        min.add(member);
	        if(min.radius<minDist){
	        	min.farthest = member;
	        	min.radius = minDist;
	        }
		}
		for (int i = 0; i < clusters2.size(); i++) {
			newRadius[i] = clusters2.get(i).radius;
		}
		return newRadius;
	}
	  
  // Main method. 
  // This program will  create
  // a jpeg file showing the expression data of each cluster. 
  public static void main(String[] astrArgs) throws Exception
  {
	  int N = 12;
	  KMeans kMeans = new KMeans();


	  	  kMeans.dist = new UserRequest.euclideanDistance();
	  	  List<RankedList> members =  readPSCFile( "data/human_cancer.pcl" );
	  	  List<Cluster> clusters = kMeans.initializeClusters(2, members);
		  kMeans.performClustering(clusters, members);

//		  String allList = "data/2157all.csv";
//		  String subLists = "data/isp.csv";
//		  kMeans.dist = new UserRequest.KendallDistance();
//		  Pair<RankedList,RankedList> centroids = readData(N, allList, subLists);
//		  List<RankedList> members =  new ArrayList<>();
//		  members.addAll((Collection<? extends RankedList>) Item.ranksListMap.values());
//	  List<Cluster> clusters = new ArrayList<>();
//	  Cluster masterCluster = new Cluster("masterCluster" );
//	  masterCluster.add(centroids.first);
//	  clusters.add(masterCluster);
//
//	  Cluster farthestCluster = new Cluster("farthestCluster" );
//	  farthestCluster.add(centroids.second);
//	  clusters.add(farthestCluster);
//
//	  kMeans.performClustering(clusters, members);
//	  List<Cluster> brokenUp = breakUpSpecifiedCluster(masterCluster);
//	  clusters.remove(masterCluster);
//	  clusters.addAll(brokenUp);
//	  kMeans.performClustering(clusters, members);

	  System.out.println("*** cluster separations "  + " ***");
	  for (int i = 0; i < clusters.size(); i++) {
		  Cluster cluster  = clusters.get(i);
		  for (int j = i+1; j < clusters.size(); j++) {
			  System.out.println(i + " " + j + " " + kMeans.dist.distance(cluster.getCentroid(), 
					  clusters.get(j).getCentroid()));
		  }
		  cluster.createJPG( cluster.name );

	  }
	  System.out.println("*** final "  + " ***");
	  for (Cluster cluster : clusters) {
		  System.out.println(cluster);
	  }
//	  for(RankedList isp : clusters.get(1).getMembers()){
//		  System.out.println(isp);
//	  }
  }

	public static Pair<RankedList,RankedList> readData(int N, String allList, String subLists)
			throws FileNotFoundException, IOException {
		RankedList master = UserRequest.readData(N, allList, subLists);
		
		Pair<Double, Double> mean = ItemSetRankedList.calculateAndPrintMeanKendall(master, false);
		System.out.println(" Mean Kendall Hausdorf dist, " +  mean.first 
		  + " Mean DCG, " +  mean.second  + " subnets " + Item.ranksListMap.size());
		Set<String> subnets = Item.ranksListMap.keySet();
		RankedList farthest = null;
		double maxKendall = 0;
		for (String string : subnets) {
			RankedList list = (RankedList) Item.ranksListMap.get(string);
			double kendall = ((double)Item.KendallHausdorf(master,list))/Item.kendallMax;
			if(maxKendall < kendall){
				maxKendall = kendall;
				farthest = list;
			}
		}
		return new Pair<RankedList,RankedList> (master, farthest);
	}
	
private static void attemptToBreakUP(KMeans kMeans, List<Cluster> clusters1, List<RankedList> members ) {
	kMeans.performClustering(clusters1, members);
	 double maxSeen = kMeans.tryToBreakUpClusters(clusters1, 23.8, members);
	 //loop doesnt work, 23.88 cannot be broken up
//	 while(maxSeen>23.8 && clusters1.size()<8){
		 maxSeen = kMeans.tryToBreakUpClusters(clusters1, 23.8, members);
//	 }
}

}
