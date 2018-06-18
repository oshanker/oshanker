/**
 * 
 */
package math;

import static org.junit.Assert.*;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.math.BigDecimal;
import java.util.Arrays;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import riemann.Gram;
import riemann.Riemann;
import riemann.Riemann.GramInfo;
import riemann.Rosser;

/**
 * @author oshanker
 *
 */
public class GSeriesTest {
    final static double[][] gramE12 = new double[][] { 
        {243.77756012466054, 7551.727850863262},
        {7551.748967244364, 14859.592051701966},
        {14859.720372293501, 22167.406308059784},
        {22167.69177527207, 29475.511171552676},
        {29475.66317618007, 36783.50533179244},
        {36783.6345750175, 44091.59099492764},
        {44091.60597178437, 51399.45918516771},
        {51399.57736648067, 58707.39943153004},
        {58707.5487591064, 66015.36710472572},
        {66015.52014966156, 73323.37514290065},
        {73323.49153814616, 80631.16991578533},
        {80631.4629245602, 87938.90285891743},
        {87939.43430890366, 95247.28904843173},
        {95247.40569117655, 102555.26217338518},
        {102555.37707137888, 109863.22425382912},
        {109863.34844951064, 117171.2193838182},
        {117171.31982557183, 124479.11858461898},
        {124479.29119956246, 131787.26254932594},
        {131787.2625714825, 139094.915172985},
        {139095.233941332, 146403.03320424244},
        {146403.20530911093, 153710.8534040047},
        {153711.1766748193, 161019.06253969827},
        {161019.1480384571, 168326.92904703945},
        {168327.11940002433, 175635.03465266858},
        {175635.09075952097, 182942.91104938296},
        {182943.06211694705, 190250.94773079225},
        {190251.03347230257, 197558.89402020318},
        {197559.00482558753, 204866.80626038337},
        {204866.97617680192, 212174.8229924622},
        {212174.94752594575, 219482.70295244805},
        {219482.918873019, 226790.86141050386},
        {226790.89021802167, 234098.80158802238},
        {234098.86156095378, 241406.80651327036},
        {241406.83290181533, 248714.5783046993},
        {248714.80424060632, 256022.741415282},
    };

	/**
	 * @throws java.lang.Exception
	 */
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
        Gram.initLogVals((int)(398942/2.2));
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
	
    @Test @Ignore
    public void testX() throws Exception{
        int k0 = 1, k1=398942;
        int R = 30040;
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        double begin = Gram.gram(offset, 243.77756012466054 );
        double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
        String zerosFile = "/Users/shankero/Documents/tmp/1e12.zeros.1001_10001002";
        BufferedReader zeroIn = new BufferedReader(new FileReader(zerosFile ));
        double upperLimit = begin + (R-40)*incr;
        Rosser.ZeroInfo zeroInput = Rosser.readZeros(upperLimit, null, zeroIn, null);
        System.out.println("{" +begin + ", " + zeroInput.lastZero+ "},");
        for (int i = 0; i < 34; i++) {
            begin = upperLimit;
            incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
            upperLimit = begin + (R-40)*incr;
            zeroInput = Rosser.readZeros(upperLimit, null, zeroIn, zeroInput.zeroInput);
            System.out.println( "{" + begin + ", " + zeroInput.lastZero + "},");
       }
    }
    
	/**
	 * Test method for {@link math.GSeries#gSeries(double)}.
	 */
	@Test
	public void test1E12() throws Exception{
		int k0 = 1, k1=398942;
		//this reduces time to 14131ms from 14326ms
		//need to tune optimum value
		DataOutputStream out = null;
        File file = new File("out/gSeriesE12.dat");
        boolean output = false;
        if (output) {
              if (!file.exists()) {
                  try {
                      file.createNewFile();
                  } catch (IOException e) {
                      e.printStackTrace();
                  }
              }
              try {
                  OutputStream os = new FileOutputStream(file);
                  BufferedOutputStream bos = new BufferedOutputStream(os);
                  // create data output stream
                  out = new DataOutputStream(bos);
              } catch (FileNotFoundException e) {
                  e.printStackTrace();
              }
        }
		int R = 30040;
		long init= System.currentTimeMillis();
		BigDecimal offset = BigDecimal.valueOf(1.0E12);
		double begin = Gram.gram(offset, gramE12[34][0] );
		double incr  = 2*Math.PI/(Math.log((offset.doubleValue()+begin)/(2*Math.PI)));
		final int initialPadding = 20;
        System.out.println(incr);
        long n0 = Gram.theta(offset.add(BigDecimal.valueOf(begin), Gram.mc), Gram.mc)
                .divide(Gram.pi, Gram.mc).longValue() - 1 - initialPadding;
        System.out.println(Gram.theta(offset.add(BigDecimal.valueOf(begin), Gram.mc), Gram.mc)
                .divide(Gram.pi, Gram.mc));
		begin -= initialPadding*incr;
		GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, R);
		long end = System.currentTimeMillis();
		System.out.println("evaluateWithOffset calc for " + R + ": " + (end - init) + "ms");
        if (output) {
            out.writeLong(n0);
    		for (int i = 0; i < gAtBeta.gAtBeta.length; i++) {
                out.writeDouble(gAtBeta.gAtBeta[i][0]);
                out.writeDouble(gAtBeta.gAtBeta[i][1]);
    		}
    		out.close();
            InputStream is = new FileInputStream(file);
            
            // create buffered input stream.
            BufferedInputStream bis = new BufferedInputStream(is);

            // create data input stream to read data in form of primitives.
            DataInputStream in = new DataInputStream(bis);
            int i = 0;
            System.out.println(in.readLong());
            while (in.available() > 0) {
                double g0 = in.readDouble();
                double g1 = in.readDouble();
                 
                if(i >= initialPadding){System.out.println(  g0 + ", " +   g1);}
                if(++i > initialPadding+1){break;}
            }
            in.close();
        }
		System.out.println(gAtBeta.riemannZeta(gAtBeta.gAtBeta[initialPadding], begin));
		//g  : [-0.33143958775035764, 0.0733285174786178] 1287.5146091794
		double zero = gramE12[34][1];
		double[] gFromBLFI = gAtBeta.diagnosticBLFISumWithOffset( zero, 4, initialPadding, 1.6E-7);
		double zeta = gAtBeta.riemannZeta(gFromBLFI, zero);
		System.out.println("g  : " + Arrays.toString(gFromBLFI) + " zeta " + zeta);
		assertTrue(Math.abs(zeta) < 0.000002);

	}


	/**
	 * Compare zeta from F-series with zeta from Riemann evaluation.
     * Z(t) = Real(exp(−i*theta(t))F(1,floor(tau); t)) + R(t)
     * dtheta/dt = ln(t/(2pi))/2;
	 */
	@Test
	public void testZeroLargeOffset() {
		double[][] fAtBeta = null;
		double[] begin = {243.8749480149, 1436233.106281030331450810};
		int k0 = 1, k1=0;
		int R = 2;
		double lnsqrtArg1 = 0;
		double basetheta = 0;
		double dsqrtArg1 = 0;
		double tbase = 0;
		double basesqrtArg1 = 0;
		long offset = (long) 1.0E12;
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
			assertTrue(i + " ", Math.abs(zeta)  < 5.0E-7);
			assertTrue(i + " ", Math.abs(zeta-zetaFromRiemann)  < 5.0E-7);
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
