/**
 * 
 */
package rank;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import static rank.DiscountedCumulativeGain.dCG;

/**
 * @author oshanker
 *
 */
public class DiscountedCumulativeGainTest {

	final static double tolerance = 0.0001d;
	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@Before
	public void setUp() throws Exception {
	}

	/**
	 * @throws java.lang.Exception
	 */
	@After
	public void tearDown() throws Exception {
	}

	@Test @Ignore
	public void testRankedList() {
		List<String> masterNames = Arrays.asList("nya", "dca", "bf1", "cha", "mib", "ne1", "daa", "dnb", "gq1", "laa", "sja", "lob");
		ItemSetRankedList master = new ItemSetRankedList(masterNames.size());
		master.add(masterNames);
		System.out.println(master);
		//assertTrue(dCG(RankedList.set, master)-10.2719<tolerance);
		System.out.println(dCG(Item.set, master));
		ItemSetRankedList reversed = new ItemSetRankedList(masterNames.size());
		int i = reversed.size;
		for (Item item: Item.set) {
			reversed.add(item, i--,-1);
		}
		System.out.println(reversed);
		//assertTrue(dCG(RankedList.set, reversed)-6.2514<tolerance);
		System.out.println(dCG(Item.set, reversed));
		reversed = new ItemSetRankedList(masterNames.size());
		i = reversed.size;
		for (Item item: Item.set) {
			if(i==1){continue;}
			reversed.add(item, i--,-1);
		}
		System.out.println(reversed);
		System.out.println(dCG(Item.set, reversed));
	}
	
	@Test
	public void test() {
		double[] rel = new double[]{3,2,3,0,1,2};
		double dcg = dCG(rel);
		
		assertTrue(Math.abs(dcg - 6.861126688593501)<tolerance);
		List<Double> ranked = new ArrayList<>();
		for (int i = 0; i < rel.length; i++) {
			ranked.add(rel[i]);
		}
		dcg = dCG(ranked);
		assertTrue(Math.abs(dcg - 6.861126688593501)<tolerance);
		List<Double> ordered = new ArrayList<>(ranked);
		Collections.sort(ordered);
		System.out.println(ordered);
		dcg = dCG(ordered);
		assertTrue(Math.abs(dcg - 4.721462852745936)<tolerance);
		Collections.sort(ordered, Collections.reverseOrder());
		System.out.println(ordered);
		dcg = dCG(ordered);
		assertTrue(Math.abs(dcg - 7.140995184095701)<tolerance);
	}

}
