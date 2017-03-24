/**
 * 
 */
package ml.kmeans;

import rank.RankedList;

/**
 * @author oshanker
 *
 */
public interface Distance {
	double distance(RankedList first, RankedList other);

}
