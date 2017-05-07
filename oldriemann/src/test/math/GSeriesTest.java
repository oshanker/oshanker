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
import riemann.Riemann.GramInfo;

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
		double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
		GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, R);
		long end = System.currentTimeMillis();
		System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
		System.out.println(gAtBeta.riemannZeta(gAtBeta.gAtBeta[0], begin));
		//g  : [-0.33143958775035764, 0.0733285174786178] 1287.5146091794
		double[] gFromBLFI = gAtBeta.blfiSumWithOffset( 1287.5146091794, 4);
		double zeta = gAtBeta.riemannZeta(gFromBLFI, 1287.5146091794);
		System.out.println("g  : " + Arrays.toString(gFromBLFI) + " zeta " + zeta);
		assertTrue(Math.abs(gFromBLFI[0] - (-0.33143958775035764)) + Math.abs(gFromBLFI[1] - 0.0733285174786178) < 0.000005);

	}

	/**
	 * Test method for {@link math.GSeries#gSeries(double)}.
	 */
	@Test
	public void test1E12() {
		int k0 = 1, k1=398942;
		Gram.initLogVals(k1);
		int R = 10000;
		long init= System.currentTimeMillis();
		BigDecimal offset = BigDecimal.valueOf(1.0E12);
		double begin = 243.77756012466052947405878015472510;
		double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
		GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, R);
		long end = System.currentTimeMillis();
		System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
		System.out.println(gAtBeta.riemannZeta(gAtBeta.gAtBeta[0], begin));
		//g  : [-0.33143958775035764, 0.0733285174786178] 1287.5146091794
		double zero = 365.657441214932156903;
		double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, 20, 1.6E-7);
		double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
		System.out.println("g  : " + Arrays.toString(gFromBLFI) + " zeta " + zeta);
		assertTrue(Math.abs(zeta) < 0.000005);

	}


	/**
	 * Compare zeta from F-series with zeta from Riemann evaluation.
     * Z(t) = Real(exp(−i*theta(t))F(1,floor(tau); t)) + R(t)
     * dtheta/dt = ln(t/(2pi))/2;
	 */
	@Test
	public void testZeroLargeOffset() {
		double[][] fAtBeta = null;
		double[] begin = {192.309350419702134727, 1921602.793342093316678265};
		int k0 = 1, k1=0;
		int R = 2;
		double lnsqrtArg1 = 0;
		double basetheta = 0;
		double dsqrtArg1 = 0;
		double tbase = 0;
		double basesqrtArg1 = 0;
		long offset = (long) 1.0E15;
		for (int i = 0; i < begin.length; i++) {
			double tincr =  (begin[i]-begin[0]) ; 
			BigDecimal tval = new BigDecimal(begin[i], Gram.mc).add(
					BigDecimal.valueOf(offset), Gram.mc);
			double predictedSqrtArg1 = 0;
			double theta = 0;
			if(i == 0 ){
				dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tval.doubleValue()));
				tbase = tval.doubleValue();
				BigDecimal t2 = tval.divide(Gram.bdTWO);
				BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
				basesqrtArg1 = sqrtArg1.doubleValue();
				k1 = (int)basesqrtArg1;
				Gram.initLogVals(k1);
				BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, Gram.mc);
				lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
				fAtBeta = GSeries.fSeries(k0, k1, begin[1]-begin[0], R, tval);
				//theta should be 2.819633653651107
				theta = tval.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
						.subtract(Gram.pi8, Gram.mc).remainder(Gram.pi_2).doubleValue();
				basetheta = theta;
				predictedSqrtArg1 = basesqrtArg1 ;
				//0.24359904690398668
				GramInfo[] gram = Gram.RZGram(3945951431271L, 2, Gram.mc);
				System.out.println("gram " + Arrays.toString(gram));
			} else {
				theta = (basetheta + lnsqrtArg1*tincr
						+tincr*tincr/(4*tbase))%(2*Math.PI);
				predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
				
				BigDecimal tBase = new BigDecimal(begin[i], Gram.mc).add(
						BigDecimal.valueOf(offset), Gram.mc);
				BigDecimal alphaBD = (Gram.log(k0).add(Gram.log(k1))).divide(Gram.bdTWO, Gram.mc);
				double argalphaBase = tBase.multiply(alphaBD, Gram.mc).remainder(Gram.pi_2).doubleValue();
				double[] g = new double[2];
				g[0] = Math.cos(argalphaBase)*fAtBeta[i][0] + Math.sin(argalphaBase)*fAtBeta[i][1];
				g[1] = Math.cos(argalphaBase)*fAtBeta[i][1] - Math.sin(argalphaBase)*fAtBeta[i][0];
				System.out.println("g  : " + Arrays.toString(g) + " argalphaBase " + argalphaBase);
				
			}
			double rotatedSum = 2*( Math.cos(theta)*fAtBeta[i][0]+Math.sin(theta)*fAtBeta[i][1]);
			double correction = GSeries.correction( predictedSqrtArg1);
			double zeta = rotatedSum + correction;
			System.out.println("f  : " + Arrays.toString(fAtBeta[i])
			   + " theta " + theta + " rotatedSum " + rotatedSum
			   + " zeta " + zeta);
			double zetaFromRiemann = Riemann.riemann(begin[i], offset);
			System.out.println("zetaFromRiemann " + zetaFromRiemann);
			System.out.println("sqrtArg1[i].doubleValue() " + predictedSqrtArg1 + " correction " + correction );
			assertTrue(i + " ", Math.abs(zeta)  < 1.0E-6);
			assertTrue(i + " ", Math.abs(zeta-zetaFromRiemann)  < 2.0E-8);
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
		double[] gFromBLFI = x.blfiSumWithOffset( t0, 2);
		double[] directEval = x.gSeries(t0);
		assertTrue(Math.abs(gFromBLFI[0] - directEval[0]) + Math.abs(gFromBLFI[1] - directEval[1]) < 0.005);
		System.out.println(t0 + " sum " + gFromBLFI[0] + ", " + gFromBLFI[1] + ": " + Arrays.toString(directEval));
	}

}
