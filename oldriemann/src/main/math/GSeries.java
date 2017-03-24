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
	final int n0, n1;
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
	private final double spacing;
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
		this.n0 = n0;
		this.n1 = n1;
		coeff = new double[k1-k0+1];
		int N = n1-n0+1;
		gAtBeta = new double[N][2];
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
		for (int i = 0; i < N; i++) {
			double t = (i+n0)*spacing;
			gAtBeta[i] = gSeries(t);
		}
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
	
	public static void main(String[] args){
		int k0 = 10, k1=100;
		int N = 12;
		int minIndex = 5;
		GSeries x = new GSeries(k0, k1,minIndex, minIndex+N-1);
		double t0 = (minIndex+N/2+0.5)*x.spacing;
		System.out.println("pi/beta " + x.spacing);
		double[] sum = x.blfiSum( t0);
		System.out.println(t0 + " sum " + Arrays.toString(sum) + ": " + Arrays.toString(x.gSeries(t0)));
	}

	public  double[] blfiSum( double t0) {
		double[] sum = new double[]{0,0};
		int N = n1-n0+1;
		for (int i = 0; i < N; i++) {
			double t = (i+n0)*spacing;
			double M = 2;
			double harg = gamma*(t0-t)/M ;
			double h = Math.pow( Math.sin(harg)/harg, M);
			double sarg = beta*(t0-t) ;
			double sin = Math.sin(sarg)/sarg;
			sum[0] += gAtBeta[i][0]*h*sin;
			sum[1] += gAtBeta[i][1]*h*sin;
			System.out.println((i+n0) + "; " + t + ": " + Arrays.toString(gAtBeta[i]) );
		}
		return sum;
	}

}
