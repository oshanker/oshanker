/**
 * 
 */
package rank;

import java.io.File;

/**
 * @author oshanker
 *
 */
public interface RankedList {

	 int rank(int index);
	 int itemCount();
	 public void add(Item i, int rank, double latency);
	double[] getValues();
	String getName();
	void setName(String key);
}
