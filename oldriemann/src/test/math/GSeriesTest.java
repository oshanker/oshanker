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
import riemann.Riemann;

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
	 * Compare zeta from F-series with zeta from Riemann evaluation.
     * Z(t) = Real(exp(−i*theta(t))F(1,floor(tau); t)) + R(t)
	 */
	@Test
	public void testZeroLargeOffset() {
		double begin = 1.8475231278;
		BigDecimal tval = new BigDecimal(begin, Gram.mc).add(
				BigDecimal.valueOf(267653395647L), Gram.mc);
		
		BigDecimal t2 = tval.divide(Gram.bdTWO);
		// 206393.703762602395481916695731615017994489594 from riemann,
		BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
		//theta should be 2.819633653651107
		double theta = tval.multiply(Gram.log(sqrtArg1, Gram.mc), Gram.mc).subtract(t2, Gram.mc)
				.subtract(Gram.pi8, Gram.mc).remainder(Gram.pi_2).doubleValue();
		int k0 = 1, k1=sqrtArg1.intValue();
		Gram.initLogVals(k1);
		System.out.println("k1 " + k1);
		int R = 1;
		double[][] fAtBeta = GSeries.fSeries(k0, k1,   0.25671765, R, tval);
		double rotatedSum = 2*( Math.cos(theta)*fAtBeta[0][0]+Math.sin(theta)*fAtBeta[0][1]);
		//0.0010100624905039076
		double p = sqrtArg1.doubleValue()-k1;
		double correction = GSeries.correction(p, sqrtArg1.doubleValue(), k1);
		//sum should be -0.0010102280025219446, cf -0.0010102280024559818
		double zeta = rotatedSum + correction;
		System.out.println("f  : " + Arrays.toString(fAtBeta[0])
		   + " theta " + theta + " rotatedSum " + rotatedSum
		   + " zeta " + zeta);
		System.out.println("correction " + correction );
		double zetaFromRiemann = Riemann.riemann(begin, 267653395647L);
		System.out.println(zetaFromRiemann);
		assertTrue(Math.abs(zeta)  < 0.00001);
		assertTrue(Math.abs(zeta-zetaFromRiemann)  < 1.0E-10);
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
