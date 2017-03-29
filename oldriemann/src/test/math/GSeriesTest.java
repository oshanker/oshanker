/**
 * 
 */
package math;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import riemann.Gram;

/**
 * @author oshanker
 *
 */
public class GSeriesTest {

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

	/**
	 * Test method for {@link math.GSeries#gSeries(double)}.
	 */
	@Test
	public void testLargeOffset() {
		int k0 = 1, k1=206393;
		Gram.initLogVals(k1);
		int R = 10000;
		long init= System.currentTimeMillis();
		BigDecimal offset = BigDecimal.valueOf(267653395647L);
		double begin = 1.87383225;
		double[][] gAtBeta = GSeries.evaluateWithOffset(k0, k1, offset,  begin,  0.25671765, R);
		long end = System.currentTimeMillis();
		System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
	}

	/**
	 * Test method for {@link math.GSeries#gSeries(double)}.
	 */
	@Test
	public void testEvaluateWithOffset() {
		int k0 = 10, k1=100;
		int N = 30;
		int minIndex = 5;
		long init= System.currentTimeMillis();
		GSeries x = new GSeries(k0, k1, minIndex, minIndex+N-1);
		long end = System.currentTimeMillis();
		System.out.println("calc for " + N + ": " + (end - init) + "ms");
		BigDecimal offset = BigDecimal.valueOf(5);
		double begin = minIndex*x.spacing - offset.doubleValue();
		init= System.currentTimeMillis();
		double[][] gAtBeta = GSeries.evaluateWithOffset(k0, k1, offset,  begin,  x.spacing, N);
		end = System.currentTimeMillis();
		System.out.println("evaluateWithOffset calc for " + N + ": " + (end - init) + "ms");
		for (int i = 0; i < gAtBeta.length; i++) {
			assertTrue(Math.abs(gAtBeta[i][0] - x.gAtBeta[i][0]) + Math.abs(gAtBeta[i][1] - x.gAtBeta[i][1]) < 0.00001);
		}
	}

	/**
	 * Test method for {@link math.GSeries#gSeries(double)}.
	 */
	@Test
	public void testGSeries() {
		int k0 = 10, k1=100;
		int N = 30;
		long init= System.currentTimeMillis();
		GSeries x = new GSeries(k0, k1, 5, 5+N-1);
		long end = System.currentTimeMillis();
		System.out.println("calc for " + N + ": " + (end - init) + "ms");
		int minIndex = 5;
		double t0 = (minIndex+N/2+0.5)*x.spacing;
		double[] gFromBLFI = x.blfiSum( t0, 2);
		double[] directEval = x.gSeries(t0);
		assertTrue(Math.abs(gFromBLFI[0] - directEval[0]) + Math.abs(gFromBLFI[1] - directEval[1]) < 0.005);
		System.out.println(t0 + " sum " + gFromBLFI[0] + ", " + gFromBLFI[1] + ": " + Arrays.toString(directEval));
	}

}
