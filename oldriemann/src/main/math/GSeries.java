package math;

import java.util.Arrays;

/**
 * Calculates the G-series useful in the study of Riemann zeta
 * function.
 * zeta(1/2 + it) = exp(−i*theta(t))Z(t)
 * tau = sqrt(t/(2*pi))
 * Z(t) = Real(exp(−i*theta(t))F(1,floor(tau); t)) + R(t)
 * 
 * @author oshanker
 *
 */
public class GSeries {
	final int k0, k1;
	final int r0, r1;
	/**
	 * The coefficients in the g series.
	 * stores from k0 to k1
	 */
	final double[] coeff;
	/**
	 * The ln term in the g series.
	 * stores from k0 to k1
	 */
	final double[] ln;
	/**
	 * rotation from F to G:
	 * G(t) = exp(−i*alpha*t)F(t)
	 */
	final double alpha;
	
	/**
	 * Store g at n*beta, for n from n0 to n1 (k is k0 to k1)
	 */
	final double[][] gAtBeta;
	private final double beta;
	final double spacing;
	private final double gamma;

	/**
	 * initialize the g-series for the given range of terms.
	 * precalculates the 1/sqrt(k) and ln(k) terms.
	 * @param k0
	 * @param k1
	 */
	public GSeries(int k0, int k1, int n0, int n1) {
		this.k0 = k0;
		this.k1 = k1;
		this.r0 = n0;
		this.r1 = n1;
		coeff = new double[k1-k0+1];
		int R = n1-n0+1;
		gAtBeta = new double[R][2];
		ln = new double[k1-k0+1];
		for (int i = k0; i <= k1; i++) {
			coeff[i-k0] = 1.0/Math.sqrt(i);
			ln[i-k0] = Math.log(i);
		}
		alpha = (Math.log(k0)+ Math.log(k1))/2.0;
		double tau = Math.log(k1/k0)/2.0;
		double lambda = 2.0d;
		beta = lambda*tau;
		spacing = Math.PI/beta;
		gamma = beta -tau;
		for (int i = 0; i < R; i++) {
			double t = (i+n0)*spacing;
			gAtBeta[i] = gSeries(t);
		}
	}
	
	/**
	 * Evaluate gSeries without calling cos and sin functions more often than
	 * necessary.
	 * Offset is placeholder, it will be extended to BigDouble for large height
	 * calculations.
	 * @param k0
	 * @param k1
	 * @param offset
	 * @param begin
	 * @param incr
	 * @param R
	 * @return
	 */
	public static double[][] evaluateWithOffset(int k0, int k1, double offset, double begin, double incr, int R){
		double[][] gAtBeta = new double[R][2];
		double tBase = offset + begin;
		for (int i = k0; i <= k1; i++) {
			//evaluate one term in the series, for all t.
			double coeff = 1/Math.sqrt(i);
			double costlni = Math.cos(tBase*Math.log(i));
			double sintlni = Math.sin(tBase*Math.log(i));
			double cosdlni = Math.cos(incr*Math.log(i));
			double sindlni = Math.sin(incr*Math.log(i));
			for (int j = 0; j < R; j++) {
				gAtBeta[j][0] += coeff*costlni;
				gAtBeta[j][1] += coeff*sintlni;
				//now set values for next t
				double tmpCos = costlni*cosdlni - sintlni*sindlni;
				sintlni = sintlni*cosdlni + costlni*sindlni;
				costlni = tmpCos;
			}
		}
		double alpha = (Math.log(k0)+ Math.log(k1))/2.0;
		double costalpha = Math.cos(tBase*alpha);
		double sintalpha = Math.sin(tBase*alpha);
		double cosdalpha = Math.cos(incr*alpha);
		double sindalpha = Math.sin(incr*alpha);
		for (int j = 0; j < R; j++) {
			double tmp = gAtBeta[j][0]*costalpha + gAtBeta[j][1]*sintalpha;
			gAtBeta[j][1] = -gAtBeta[j][0]*sintalpha + gAtBeta[j][1]*costalpha;
			gAtBeta[j][0] = tmp;
			//now set values for next t
			double tmpCos = costalpha*cosdalpha - sintalpha*sindalpha;
			sintalpha = sintalpha*cosdalpha + costalpha*sindalpha;
			costalpha = tmpCos;
		}
		
		return gAtBeta;
	}
	
	/**
	 * Calculate the gSeries for arg t.
	 * @param t
	 * @return
	 */
	double[] gSeries(double t){
		double[] g = new double[2];
		double f0 = 0, f1 = 0;
		for (int i = k0; i <= k1; i++) {
			f0 += coeff[i-k0]*Math.cos(t*ln[i-k0]);
			f1 += coeff[i-k0]*Math.sin(t*ln[i-k0]);
		}
		g[0] = Math.cos(alpha*t)*f0 + Math.sin(alpha*t)*f1;
		g[1] = Math.cos(alpha*t)*f1 - Math.sin(alpha*t)*f0;
		return g;
	}

	public  double[] testblfiSum( double t0, int M) {
		double[] directEval = gSeries(t0);
		double[] sum = new double[]{0,0};
		int midIdx = (int) (t0/spacing);
		SUM:
		for (int term = 0; term < 100; term++) {
			for (int j = 0; j < 2; j++) {
				int i = midIdx + (j==0?-term:(term+1));
				if(i>r1 || i < r0){
					System.out.println("breaking " + i);
					break SUM;}
				double t = (i)*spacing;
				double harg = gamma*(t0-t)/M ;
				double h = Math.pow( Math.sin(harg)/harg, M);
				double sarg = beta*(t0-t) ;
				double sin = Math.sin(sarg)/sarg;
				sum[0] += gAtBeta[i-r0][0]*h*sin;
				sum[1] += gAtBeta[i-r0][1]*h*sin;
			}
			System.out.println((term) + " : " + Arrays.toString(sum) 
			   + " : " + (Math.abs(sum[0] - directEval[0]) + " : " + Math.abs(sum[1] - directEval[1])));
		}
		return sum;
	}

	public  double[] blfiSum( double t0, int M) {
		double[] sum = new double[]{0,0};
		int midIdx = (int) (t0/spacing);
		SUM:
		for (int term = 0; term < 8; term++) {
			for (int j = 0; j < 2; j++) {
				int i = midIdx + (j==0?-term:(term+1));
				if(i>r1 || i < r0){break SUM;}
				double t = (i)*spacing;
				double harg = gamma*(t0-t)/M ;
				double h = Math.pow( Math.sin(harg)/harg, M);
				double sarg = beta*(t0-t) ;
				double sin = Math.sin(sarg)/sarg;
				sum[0] += gAtBeta[i-r0][0]*h*sin;
				sum[1] += gAtBeta[i-r0][1]*h*sin;
			}
		}
		return sum;
	}
	
	public static void main(String[] args){
		int k0 = 10, k1=100;
		int N = 25;
		int minIndex = 5;
		GSeries x = new GSeries(k0, k1,minIndex, minIndex+N-1);
		System.out.println("pi/beta " + x.spacing);
		double[] offsets = { 0.3, 0.5, 0.7};
		for (int i = 0; i < offsets.length; i++) {
			double t0 = (minIndex+N/2+offsets[i])*x.spacing;
			double[] sum = x.testblfiSum( t0, 3);
			System.out.println(t0 + " sum " + Arrays.toString(sum) + ": " + Arrays.toString(x.gSeries(t0)));
		}
	}

}
