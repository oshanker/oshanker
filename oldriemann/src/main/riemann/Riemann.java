package riemann;

import java.math.BigDecimal;
import java.math.MathContext;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;


/*
 Author: O. Shanker.
 oshanker At gmail dot com
 http://sites.google.com/site/oshanker/
 
 Brute force evaluation of Riemann function, using Riemann Siegsl series.
 Useful for checking output of more optimized approaches.


 # Permission is hereby granted, free of charge, to any person obtaining 
 # a copy of this software and associated documentation files (the
 # "Software"), to deal in the Software without restriction, including
 # without limitation the rights to use, copy, modify, merge, publish,
 # distribute copies of the Software, and to
 # permit persons to whom the Software is furnished to do so, subject to
 # the following conditions:
 # 
 # The above copyright notice and this permission notice shall be
 # included in all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 # EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 # MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 # IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 # CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 # TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 # SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
public class Riemann {

	static boolean verbose = true;
	static double[] baseArgI = null;
	static int oldN = 0;
	static BigDecimal baseT = null;

	final static Map<String, BigDecimal> offsets = new ConcurrentHashMap<String, BigDecimal>();

	static double[] c1coeff = { -.02682510262837534703, +.01378477342635185305,
			.03849125048223508223, +.00987106629906207647,
			-.00331075976085840433, -.00146478085779541508,
			-.00001320794062487696,
	};
	
	static NumberFormat nf = NumberFormat.getInstance();
	static {
		nf.setMinimumFractionDigits(12);
		nf.setMaximumFractionDigits(12);
		nf.setGroupingUsed(false);
		offsets.put("12", BigDecimal.valueOf(  267653395647L));// 12
		offsets.put("A", BigDecimal.valueOf(  7954022502370L));// A
		offsets.put("B", BigDecimal.valueOf(   323393653040L));// B
		offsets.put("C", BigDecimal.valueOf( 18580341990000L));// C
		offsets.put("D", BigDecimal.valueOf( 18523741991600L));//
		offsets.put("E", BigDecimal.valueOf(  4307762397300L));//
		offsets.put("F", BigDecimal.valueOf( 16671318581100L));//
		offsets.put("G", BigDecimal.valueOf(  8847150598000L));//
		offsets.put("H", BigDecimal.valueOf( 10123949469500L));//
		offsets.put("I", BigDecimal.valueOf( 20825125156950L));//
		offsets.put("J", BigDecimal.valueOf( 17335277232200L));//
		offsets.put("K", BigDecimal.valueOf(  5907264585800L));//
		offsets.put("L", BigDecimal.valueOf( 21285800773550L));//
		offsets.put("13", BigDecimal.valueOf( 1034741742900L));// 13
		offsets.put("M", BigDecimal.valueOf(  2124447368580L));// M
		offsets.put("N", BigDecimal.valueOf(  2414113624150L));
		offsets.put("O", BigDecimal.valueOf(    10854395950L));
		offsets.put("P", BigDecimal.valueOf(   694884271750L));
		offsets.put("Q", BigDecimal.valueOf(  1336685304900L));
		offsets.put("W", BigDecimal.valueOf(    35615956500L));
		offsets.put("Lune", BigDecimal.valueOf(   388858870L));
		offsets.put("Z4b", BigDecimal.valueOf(1674556748200L));
		offsets.put("Z5", BigDecimal.valueOf(  935203331160L));
	}

	private static void initN(MathContext mc, BigDecimal t, int N) {
		if(N == oldN) { return; }
		oldN = N;
		Gram.initLogVals(N);
		baseArgI = new double[N];
		baseT = t;
		for (int i = 1; i <= N; i++) {
			baseArgI[i - 1] = (t.multiply(Gram.log(i), mc)).remainder(Gram.pi_2)
					.doubleValue();
		}
	}

	/**
	 * The correction term in the Riemann Siegel series.
	 * @param mc
	 * @param sqrtArg1
	 * @param N
	 * @return
	 */
	public static double rTerm(MathContext mc, BigDecimal sqrtArg1, int N) {
		double p = sqrtArg1.remainder(BigDecimal.ONE).doubleValue();
		double fourthRoot = Gram.sqrt(sqrtArg1, mc, 1.0E-18).doubleValue();
		double c0formula = Math.cos(2 * Math.PI * (p * p - p - 1.0 / 16))
				/ (fourthRoot * Math.cos(2 * Math.PI * p));
		double c1 = c1coeff(fourthRoot, p);
		double R = c0formula + c1;
		if (N % 2 == 0) {
			R = -R;
		}
		return R;
	}

	private static void bracketZeros(BigDecimal offset) {
		double[] vals = new double[] { 
// odlyzko 12
				//			     1.8475231278,
//			     2.3623669687,
//			     2.6816309165
				
			       244.02115917156451846579279784626171,
			       244.26475821846850509076039309471551,
			       244.50835726537248948123327538095505,
			       244.75195631227647157107558996548808,
				// http://sage.math.washington.edu/home/hiaryg/page/index.html 12
			       244.158906912980683962,
			       244.367502584863394599,
//			       244.588579452072075626,
//			       244.920599505825861697,
//			       245.080792332878812378,
//			       245.370714305188583892
				
				// http://sage.math.washington.edu/home/hiaryg/page/index.html 15
	             //  192.309350419702134727,
				 };
		System.out.println("bracketZeros " + offset);
		BigDecimal[] tvals = new BigDecimal[vals.length];
		for (int i = 0; i < vals.length; i++) {
			tvals[i] = offset.add(BigDecimal.valueOf(vals[i]), Gram.mc);
		}
		long idx = 0;
		for (int i = 0; i < 1; i++) {
			BigDecimal heck = Gram.theta(tvals[i], Gram.mc)
					.divide(Gram.pi, Gram.mc);
			idx = heck.longValue();
			System.out.println("val " + tvals[i] + " theta/pi " + heck + ", lower gram idx " + idx);
		}
		GramInfo[] x = Gram.RZGram(idx, 4, Gram.mc);
		for (GramInfo gramInfo : x) {
			System.out.println("val-offset " + gramInfo.grampt.subtract(offset) + " idx " + gramInfo.idx);
		}
		System.out.println("diff " + x[x.length-1].grampt.subtract(x[0].grampt).doubleValue()/(x.length-1));
		long init = System.currentTimeMillis();
		double[] riemann = riemann(tvals, Gram.mc);
		for (int i = 0; i < vals.length; i++) {
			System.out.println(tvals[i].subtract(offset, Gram.mc) + "  Z "
					+ riemann[i]);
		}
		long end = System.currentTimeMillis();
		System.out.println("calc for " + riemann.length + " " + (end - init) + "ms");
	}

	/**
	 * The c1 correction term in the Riemann Siegel series.
	 * @param fourthRoot
	 * @param p
	 * @return
	 */
	public static double c1coeff(double fourthRoot, double p) {
		double mult = (2 * p - 1);
		double term = mult;
		double c1 = 0;
		mult = mult * mult;
		for (int i = 0; i < c1coeff.length; i++) {
			c1 += c1coeff[i] * term;
			term *= mult;
		}
		c1 /= Math.pow(fourthRoot, 3);
		return c1;
	}

	static BigDecimal myBigTheta(BigDecimal t, MathContext mc) {
		BigDecimal arg1 = Gram.sqrt(t.divide(Gram.pie2, mc), mc, 1.0E-15);
		BigDecimal arg2 = Gram.sqrt(arg1, mc, 1.0E-15);
		Gram.initLogVals(arg2.intValue());
		BigDecimal theta = t.multiply(Gram.log(arg2,mc),mc).multiply(Gram.bdTWO)
				.subtract(Gram.pi8 ,mc);
		return theta;
	}	
	
	public static void main(String[] args) throws Exception {
//		bracketZeros(offsets.get("12"));
		bracketZeros(new BigDecimal("1000000000000"));

	}

	/** calculate zeta function for (t+offset) */
	public static double riemann(double t, long offset) {

		BigDecimal[] bdtz;
		try {
			bdtz = new BigDecimal[] { new BigDecimal(t, Gram.mc).add(
					BigDecimal.valueOf(offset), Gram.mc) };

		} catch (Exception e) {
			System.out.println("t " + t + ", " + e.getMessage());
			throw new RuntimeException(e);
		}

		return riemann(bdtz, Gram.mc)[0];
	}

	/** calculate zeta function for set of arguments. */
	public static double[] riemann(BigDecimal[] t, MathContext mc) {
		double[] riemann = new double[t.length];
		for (int j = 0; j < riemann.length; j++) {
			BigDecimal tval = t[j];
			BigDecimal t2 = tval.divide(Gram.bdTWO, mc);
			//BigDecimal arg1 = t2.divide(Gram.pi, mc);
			BigDecimal sqrtArg1 = Gram.sqrt(t2.divide(Gram.pi, mc), Gram.mc, 1.0E-21);
			double theta = tval.multiply(Gram.log(sqrtArg1, mc), mc).subtract(t2, mc)
					.subtract(Gram.pi8, mc).remainder(Gram.pi_2).doubleValue();
			int N = sqrtArg1.intValue();
			double correction = rTerm(mc, sqrtArg1, N);
			double sum = 0;
			if (oldN != N) {
				initN(mc, tval, N);
				for (int i = 1; i <= N; i++) {
					sum += Math.cos(theta - baseArgI[i - 1]) / Math.sqrt(i );
				}
			} else {
				double del = tval.subtract(baseT, mc).doubleValue();
				for (int i = 1; i <= N; i++) {
					double argi = baseArgI[i - 1] + del * Gram.logdbl(i);
					double rcos = Math.cos(theta - argi) / Math.sqrt(i );
					sum += rcos;
				}
			}
			sum = 2 * sum;

			riemann[j] = sum + correction;
		}
		return riemann;
	}

	public static double riemannTestEval(BigDecimal t, MathContext mc) {
		BigDecimal t2 = t.divide(Gram.bdTWO, mc);
		BigDecimal arg1 = t2.divide(Gram.pi, mc);
		BigDecimal sqrtArg1 = Gram.sqrt(arg1, Gram.mc, 1.0E-18);
		int N = sqrtArg1.intValue();
		double R = rTerm(mc, sqrtArg1, N);
		//initSqrt(N);
		Gram.initLogVals(N);
		BigDecimal term = BigDecimal.ONE.divide(
				t.multiply(BigDecimal.valueOf(48), mc), mc);
		BigDecimal theta = t.multiply(Gram.log(sqrtArg1, mc), mc).subtract(t2, mc)
				.subtract(Gram.pi8, mc).add(term, mc);
		double sum = 0;
		for (int i = 1; i <= N; i++) {
			double arg = theta.subtract(t.multiply(Gram.log(i), mc), mc)
					.remainder(Gram.pi_2).doubleValue();
			sum += Math.cos(arg) / Math.sqrt(i );
		}
		sum = 2 * sum;

		double riemann = sum + R;
		return riemann;
	}

	/**
	 * Stores a t and the z value for the t.
	 * @author oshanker
	 *
	 */
	public static class RZpoint implements Comparable<RZpoint> {
		final double t;
		double z;
		BigDecimal bdOffset;

		public RZpoint(double t) {
			this.t = t;
		}

		public RZpoint(double t, double z) {
			this.t = t;
			this.z = z;
		}

		public RZpoint(GramInfo gramInfo, BigDecimal bdOffset) {
			this.t = gramInfo.grampt.subtract(bdOffset, Gram.mc).doubleValue();
			this.z = gramInfo.zeta;
			this.bdOffset = bdOffset;
		}

		public String toString() {
			return "RZpoint [t=" + t + ", z=" + z + "]";
		}

		public int compareTo(RZpoint arg0) {
			RZpoint other = (RZpoint) arg0;
			if (t < other.t) {
				return -1;
			}
			if (t == other.t) {
				return 0;
			}
			return 1;
		}

		public int hashCode() {
			final int prime = 31;
			int result = 1;
			long temp;
			temp = Double.doubleToLongBits(t);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(z);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			return result;
		}

		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			RZpoint other = (RZpoint) obj;
			if (Double.doubleToLongBits(t) != Double.doubleToLongBits(other.t))
				return false;
			if (Double.doubleToLongBits(z) != Double.doubleToLongBits(other.z))
				return false;
			return true;
		}
	}

	public static class GramInfo implements Comparable<GramInfo>{
		public final BigDecimal grampt;
		double zeta;
		double[] cos = new double[10];
		double[] sin = new double[9];
		double der;
		public final BigDecimal idx;
		RZpoint point;

		public GramInfo(BigDecimal grampt, long idx) {
			this.grampt = grampt;
			this.idx = BigDecimal.valueOf(idx);
		}

		public GramInfo(BigDecimal grampt, BigDecimal idx) {
			this.grampt = grampt;
			this.idx = idx;
		}

		RZpoint setRZpoint(BigDecimal bdOffset) {
			if (point == null) {
				point = new RZpoint(this, bdOffset);
			}
			return point;
		}

		static TreeSet<RZpoint> toRZPoint(Iterable<GramInfo> riemann,
				BigDecimal bdOffset) {
			TreeSet<RZpoint> gram = new TreeSet<RZpoint>();
			for (GramInfo gramInfo : riemann) {
				gram.add(gramInfo.setRZpoint(bdOffset));
			}
			return gram;
		}

		static GramInfo fromStringArray(String[] parsed, BigDecimal bdOffset) {
			BigDecimal grampt = BigDecimal.valueOf(
					Double.parseDouble(parsed[0])).add(bdOffset, Gram.mc);
			GramInfo gramInfo = new GramInfo(grampt, Long.parseLong(parsed[22]));
			gramInfo.zeta = Double.parseDouble(parsed[1]);
			gramInfo.der = Double.parseDouble(parsed[21]);
			int offset = 2;
			for (int j = 0; j < gramInfo.cos.length; j++) {
				gramInfo.cos[j] = Double.parseDouble(parsed[j + offset]);
			}
			offset += gramInfo.cos.length;
			for (int j = 0; j < gramInfo.sin.length; j++) {
				gramInfo.sin[j] = Double.parseDouble(parsed[j + offset]);
			}
			return gramInfo;
		}

		public String toString(BigDecimal offset) {
			String ret = nf.format(offset==null?grampt:grampt.subtract(offset)) + " \t" + zeta
					+ " \t" + Arrays.toString(cos) + " \t"
					+ Arrays.toString(sin) + " \t" + der + " \t" + idx;
			ret = ret.replaceAll("[\\[\\],]", "");
			return ret;
		}

		@Override
		public String toString() {
			return "GramInfo [grampt=" + grampt + ", idx=" + idx + "]";
		}

		public int compareTo(GramInfo arg0) {
			GramInfo other = (GramInfo) arg0;
			return idx.compareTo(other.idx);
		}
	}
}
