// This class represents a cluster of subnets with similar latency data. It
// contains methods for creating an image of the cluster's subnets data and
// for returning the centroid of the cluster.
package ml.kmeans;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;
import javax.imageio.ImageIO;

import rank.RankedList;

import java.io.File;

public class Cluster 
{
  // Constants for creating Jpeg images. Do not change these.
  public static double C_D_CUTOFF = 2;
  private static final int C_IX  = 8;
  private static final int C_IY  = 4;
  
  // Data members
  private HashSet<RankedList> members = new HashSet<>(); // List of subnets in cluster
   public String name;
   public double radius = 0;
   public RankedList farthest;
   public UserRequest centroid;
  
  // Constructor; creates an empty cluster.
  public Cluster( String id) 
  {
	  this.name = id;
  }
  
  public UserRequest getCentroid() {
	return centroid;
}

public HashSet<RankedList> getMembers() {
	return members;
}

public String getName() {
	return name;
}

public double getRadius() {
	return radius;
}

public RankedList getFarthest() {
	return farthest;
}

  public int size(){
	  return members.size();
  }
  
  // Prints the names of all genes in the cluster
  public void printNames() 
  {
    System.out.println(members);
  }
  
  // Adds a member to the cluster
  public void add(RankedList rankedList) 
  {
    // Append a gene to the set
    members.add(rankedList); 
  }
  
  // Adds a member to the cluster
  public void addAll(Collection<UserRequest> gene) 
  {
    // Append a gene to the set
    members.addAll(gene); 
  }
  
  public void clear()
  {
    members.clear();
    radius = 0;
    farthest = null;
  }

  public boolean remove(UserRequest member)
  {
     return members.remove(member);
  }

  // Returns the centroid of the cluster 
  public UserRequest generateCentroid() 
  {
	  int length = members.iterator().next().itemCount();
	  double[] latencies = new double[length];
	  for (RankedList userRequest : members) {
		  for (int i = 0; i < length; i++) {
			  latencies[i] += userRequest.getValues()[i];
		  }

	  }
	  for (int i = 0; i < length; i++) {
		  latencies[i] /= members.size();
	  }
	  centroid = new UserRequest(name, latencies);
      return centroid; 
  }
  
  @Override
public String toString() {
	return name + " " + size() + " " + radius;
}

//Creates an image of this cluster's expression data. The image will be
 // stored in file "<fileName><id>.jpg". Do not change this method.
 public void myCreateJPG(String fileName) 
 {
   String   strOut;
   int    i, j, iGenes, iConditions;
   double   dValue;
   BufferedImage bimg;
   Graphics2D  gr2d;
   Color   colr;
   
   strOut = "out/" + fileName  + ".jpg";
   
   // Initialize some values
   Vector<RankedList> geneVector = new Vector<RankedList>();
   Iterator<RankedList> it = members.iterator();
   while (it.hasNext())
       geneVector.add(it.next());

   Comparator<RankedList> GENE_ORDER = new Comparator<RankedList>() {
      public int compare(RankedList g1, RankedList g2) {
         return g1.getName().compareTo(g2.getName());
      }
   };

   Collections.sort(geneVector, GENE_ORDER);

   iGenes = geneVector.size();
   iConditions = (geneVector.get(0)).getValues().length;
   
   // Create the empty image and graphics2D
   bimg = new BufferedImage(C_IX * iConditions, C_IY * iGenes,
                            BufferedImage.TYPE_INT_RGB);
   gr2d = bimg.createGraphics();
   
   // Draw a rectangle for each gene/condition pair
   for (i = 0; i < iGenes; ++i) 
   {
     for (j = 0; j < iConditions; ++j) 
     {
       dValue = (geneVector.get(i)).getValues()[j];
       if (dValue < 0){
         throw new IllegalStateException("ilegal latency " + dValue);
       }
       else {
    	   float r  = 0, g = 0, b = 0;
    	   if(dValue >= 175){
    		   r = 1.0f;
    	   } else if(dValue <= 75){
    		   g = (float) ((200 -dValue)/200);
    	   } else {
    		   b = 1.0f;
    	   }
           colr = new Color(r,g,b);
       }
       gr2d.setColor(colr);
       gr2d.fill(new Rectangle2D.Float(j * C_IX, i * C_IY, C_IX, C_IY)); 
     } 
   }
   gr2d.dispose();
   
   try 
   {
     // Output the image to file
     File outFile = new File(strOut);
     ImageIO.write(bimg, "jpg", outFile);
   }
   catch (IOException e) 
   {
     System.out.println("ERROR: Unable to write image in " + strOut + "."); 
   } 
 }

// Creates an image of this cluster's expression data. The image will be
  // stored in file "<fileName><id>.jpg". Do not change this method.
  public void createJPG(String fileName) 
  {
    String   strOut;
    int    i, j, iGenes, iConditions;
    double   dValue;
    BufferedImage bimg;
    Graphics2D  gr2d;
    Color   colr;
    
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
	
    strOut = "out/" + fileName  + ".jpg";
    
    // Initialize some values
    Vector<RankedList> geneVector = new Vector<RankedList>();
    Iterator<RankedList> it = members.iterator();
    while (it.hasNext())
        geneVector.add(it.next());

    Comparator<RankedList> GENE_ORDER = new Comparator<RankedList>() {
       public int compare(RankedList g1, RankedList g2) {
          return g1.getName().compareTo(g2.getName());
       }
    };

    Collections.sort(geneVector, GENE_ORDER);

    iGenes = geneVector.size();
    iConditions = (geneVector.get(0)).getValues().length;
    
    // Create the empty image and graphics2D
    bimg = new BufferedImage(C_IX * iConditions, C_IY * iGenes,
                             BufferedImage.TYPE_INT_RGB);
    gr2d = bimg.createGraphics();
    
    // Draw a rectangle for each gene/condition pair
    for (i = 0; i < iGenes; ++i) 
    {
      for (j = 0; j < iConditions; ++j) 
      {
        dValue = (geneVector.get(i)).getValues()[j];
        if (dValue < 0)
          colr = new Color(0.0f, (dValue < -C_D_CUTOFF) ? 1.0f
                               : (float) (dValue / -C_D_CUTOFF), 0.0f);
        else
          colr = new Color((dValue > C_D_CUTOFF) ? 1.0f
                               : (float) (dValue / C_D_CUTOFF), 0.0f, 0.0f);
        gr2d.setColor(colr);
        gr2d.fill(new Rectangle2D.Float(j * C_IX, i * C_IY, C_IX, C_IY)); 
      } 
    }
    gr2d.dispose();
    
    try 
    {
      // Output the image to file
      File outFile = new File(strOut);
      ImageIO.write(bimg, "jpg", outFile);
    }
    catch (IOException e) 
    {
      System.out.println("ERROR: Unable to write image in " + strOut + "."); 
    } 
  }
}
