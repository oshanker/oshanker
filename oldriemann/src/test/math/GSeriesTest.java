/**
 * 
 */
package math;

import static org.junit.Assert.*;

import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

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
	public void testGSeries() {
		int k0 = 10, k1=100;
		GSeries x = new GSeries(k0, k1);
		double tau = Math.log(k1/k0)/2.0;
		double lambda = 2.0d;
		double beta = lambda*tau;;
		double spacing = Math.PI/beta;
		double gamma = beta -tau;
		int N = 12;
		int minIndex = 5;
		double t0 = (minIndex+N/2+0.5)*spacing;
		System.out.println("pi/beta " + spacing);
		double[] sum = new double[]{0,0};
		for (int i = 0; i < N; i++) {
			double t = (i+minIndex)*spacing;
			double[] val = x.gSeries(t);
			double M = 2;
			double harg = gamma*(t0-t)/M ;
			double h = Math.pow( Math.sin(harg)/harg, M);
			double sarg = beta*(t0-t) ;
			double sin = Math.sin(sarg)/sarg;
			sum[0] += val[0]*h*sin;
			sum[1] += val[1]*h*sin;
			System.out.println((i+minIndex) + "; " + t + ": " + Arrays.toString(val) );
		}
		double[] directEval = x.gSeries(t0);
		assertTrue(Math.abs(sum[0] - directEval[0]) + Math.abs(sum[1] - directEval[1]) < 0.005);
		System.out.println(t0 + " sum " + sum[0] + ", " + sum[1] + ": " + Arrays.toString(directEval));
	}

}
