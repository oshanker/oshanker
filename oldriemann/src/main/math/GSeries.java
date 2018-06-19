package math;

import java.math.BigDecimal;
import java.util.Arrays;

import riemann.Gram;
import riemann.Riemann;

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
	/**
	 * rotation from F to G:
	 * G(t) = exp(−i*alpha*t)F(t)
	 */
	final double alpha;
	
	/**
	 * Store g at n*beta (k is k0 to k1)
	 */
	double[][] gAtBeta;
	private final double beta;
	final double spacing;
	private final double gamma;
	private double argalphaBase;
	private final double begin;
	final double  dsqrtArg1;
	double tbase;
	double basesqrtArg1;
	double lnsqrtArg1;
	double basetheta;
	
    public GSeries(int k0, int k1, BigDecimal offset, double begin, double incr, double[][] gAtBeta){
        this(k0, k1, offset, begin, incr);
        this.gAtBeta = gAtBeta;
    }

	public GSeries(int k0, int k1, BigDecimal offset, double begin, double incr, int R){
		this(k0, k1, offset, begin, incr);
        BigDecimal tBaseBD = new BigDecimal(begin, Gram.mc).add(
                offset, Gram.mc);
        gAtBeta = evaluateWithOffset(k0, k1, begin, incr, R, tBaseBD);
	}

    private  GSeries(int k0, int k1, BigDecimal offset, double begin, double incr) {
        this.k0 = k0;
		this.k1 = k1;
		this.begin = begin;
		this.spacing = incr;
		beta = Math.PI/spacing;
		BigDecimal tBaseBD = new BigDecimal(begin, Gram.mc).add(
				offset, Gram.mc);
		BigDecimal alphaBD = (Gram.log(k0).add(Gram.log(k1))).divide(Gram.bdTWO, Gram.mc);
		alpha = alphaBD.doubleValue();
		argalphaBase = tBaseBD.multiply(alphaBD, Gram.mc).remainder(Gram.pi_2).doubleValue();
		
		/////
		dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*tBaseBD.doubleValue()));
		tbase = tBaseBD.doubleValue();
		BigDecimal t2 = tBaseBD.divide(Gram.bdTWO);
		BigDecimal sqrtArg1 = Gram.sqrt(tBaseBD.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-21);
		basesqrtArg1 = sqrtArg1.doubleValue();
		BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, Gram.mc);
		lnsqrtArg1 = lnsqrtArg1BD.doubleValue();
		basetheta = tBaseBD.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
				.subtract(Gram.pi8, Gram.mc).remainder(Gram.pi_2).doubleValue();
		double tau = (Math.log(k1) - Math.log(k0))/2.0;
		gamma = beta -tau;
    }

	/**
	 * initialize the g-series for the given range of terms.
	 * precalculates the 1/sqrt(k) and ln(k) terms.
	 * @param k0
	 * @param k1
	 */
	GSeries(int k0, int k1, int n0, int n1) {
		this.k0 = k0;
		this.k1 = k1;
		int R = n1-n0+1;
		gAtBeta = new double[R][2];
		alpha = (Math.log(k0)+ Math.log(k1))/2.0;
		double tau = (Math.log(k1) - Math.log(k0))/2.0;
		double lambda = 2.0d;
		beta = lambda*tau;
		spacing = Math.PI/beta;
		this.begin = spacing*n0;
		gamma = beta -tau;
		dsqrtArg1 = 1.0/(2*Math.sqrt(2*Math.PI*n0*spacing));
		for (int i = 0; i < n1-n0+1; i++) {
			double t = (i+n0)*spacing;
			gAtBeta[i] = gSeries(t);
		}
	}
	
	public double riemannZeta(double[] g, double tincr){
		tincr -= begin;
		double theta = (basetheta + lnsqrtArg1*tincr
				+tincr*tincr/(4*tbase))%(2*Math.PI);
		double predictedSqrtArg1 = basesqrtArg1 + dsqrtArg1*tincr;
		double[] fAtBeta = new double[2];
		double argalphat = (argalphaBase + alpha*tincr)%(2*Math.PI);
		double cos = Math.cos(argalphat);
		double sin = Math.sin(argalphat);
		//calculate f from g    
		fAtBeta[0] = cos*g[0] - sin*g[1];
		fAtBeta[1] = sin*g[0] + cos*g[1];
		double rotatedSum = 2*( Math.cos(theta)*fAtBeta [0]+Math.sin(theta)*fAtBeta[1]);
		double correction = GSeries.correction( predictedSqrtArg1);
		double zeta = rotatedSum + correction;
		return zeta;
	}
	
	public static double correction( double sqrtArg1 ) {
		long N = (long)sqrtArg1;
		double p = sqrtArg1-N;
		double fourthRoot = Math.sqrt(sqrtArg1);
		double c0formula = Math.cos(2 * Math.PI * (p * p - p - 1.0 / 16))
				/ (fourthRoot * Math.cos(2 * Math.PI * p));
		double c1 = Riemann.c1coeff(fourthRoot, p);
		double R = c0formula + c1;
		if (N % 2 == 0) {
			R = -R;
		}
		return R;
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
	private final  double[][] evaluateWithOffset(int k0, int k1, double begin, double incr, int R, BigDecimal tBase){
		double[][] gAtBeta = fSeries(k0, k1, incr, R, tBase);
		//now ratate the f to g.
		double costalpha = Math.cos(argalphaBase);
		double sintalpha = Math.sin(argalphaBase);
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

	public static double[][] fSeries(int k0, long k1, double incr, int R, BigDecimal tBase) {
		double[][] fAtBeta = new double[R][2];
		for (int i = k0; i <= k1; i++) {
			//evaluate one term in the series, for all t.
			double coeff = 1/Math.sqrt(i);;
			double argi = tBase.multiply(Gram.log(i), Gram.mc).remainder(Gram.pi_2).doubleValue();
			double costlni = Math.cos(argi);
			double sintlni = Math.sin(argi);
			//this speeds up, but do we lose accuracy?
			// no, that is not culprit: R costlni -0.9949659697160186, cf -0.9949659697160186
			double cosdlni = Math.cos(incr*Math.log(i));
			double sindlni = Math.sin(incr*Math.log(i));
			for (int j = 0; j < R; j++) {
				fAtBeta[j][0] += coeff*costlni;
				fAtBeta[j][1] += coeff*sintlni;
				//now set values for next t
				double tmpCos = costlni*cosdlni - sintlni*sindlni;
				sintlni = sintlni*cosdlni + costlni*sindlni;
				costlni = tmpCos;
			}
		}
		return fAtBeta;
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
			f0 += Math.cos(t*Math.log(i))/Math.sqrt(i);
			f1 += Math.sin(t*Math.log(i))/Math.sqrt(i);
		}
		g[0] = Math.cos(alpha*t)*f0 + Math.sin(alpha*t)*f1;
		g[1] = Math.cos(alpha*t)*f1 - Math.sin(alpha*t)*f0;
		return g;
	}

	private  double[] testblfiSum( double t0, int M) {
		double[] directEval = gSeries(t0);
		double[] sum = new double[]{0,0};
		//int midIdx = (int) (t0/spacing);
		int midIdx = (int) ((t0-begin)/spacing);
		SUM:
		for (int term = 0; term < 100; term++) {
			for (int j = 0; j < 2; j++) {
				int i = midIdx + (j==0?-term:(term+1));
				if(i>=gAtBeta.length || i < 0){
					System.out.println("breaking " + i);
					break SUM;}
				double t = begin + (i)*spacing;
				double harg = gamma*(t0-t)/M ;
				double h = Math.pow( Math.sin(harg)/harg, M);
				double sarg = beta*(t0-t) ;
				double sin = Math.sin(sarg)/sarg;
				sum[0] += gAtBeta[i][0]*h*sin;
				sum[1] += gAtBeta[i][1]*h*sin;
			}
			System.out.println((term) + " : " + Arrays.toString(sum) 
			   + " : " + (Math.abs(sum[0] - directEval[0]) + " : " + Math.abs(sum[1] - directEval[1])));
		}
		return sum;
	}

	public  double[] diagnosticBLFISumWithOffset( double t0, int M, int terms, double tolerance) {
		double[] sum = new double[]{0,0};
		double[] oldSum = new double[]{0,0};
		int midIdx = (int) ((t0-begin)/spacing);
		SUM:
		for (int term = 0; term < terms; term++) {
			oldSum[0] = sum[0];
			oldSum[1] = sum[1];
			for (int j = 0; j < 2; j++) {
				int i = midIdx + (j==0?-term:(term+1));
				if(i>=gAtBeta.length || i < 0){
					System.out.println("breaking " + i);
					break SUM;
				}
				double t = begin + (i)*spacing;
				double harg = gamma*(t0-t)/M ;
				double h = Math.pow( Math.sin(harg)/harg, M);
				double sarg = beta*(t0-t) ;
				double sin = Math.sin(sarg)/sarg;
				sum[0] += gAtBeta[i][0]*h*sin;
				sum[1] += gAtBeta[i][1]*h*sin;
			}
			double change = Math.abs(sum[0] - oldSum[0]) +  Math.abs(sum[1] - oldSum[1]);	
			if(term%5==0 ||(change < tolerance)){
				System.out.println((term) + " : " + Arrays.toString(sum) 
				   + " : " + (Math.abs(sum[0] - oldSum[0]) + " : " + Math.abs(sum[1] - oldSum[1])));
			}
			if(change < tolerance){break;}
		}
		return sum;
	}

	double[] blfiSumWithOffset( double t0, int M) {
		double[] sum = new double[]{0,0};
		int midIdx = (int) ((t0-begin)/spacing);
		SUM:
		for (int term = 0; term < 8; term++) {
			for (int j = 0; j < 2; j++) {
				int i = midIdx + (j==0?-term:(term+1));
				if(i>=gAtBeta.length || i < 0){break SUM;}
				double t = begin + (i)*spacing;
				double harg = gamma*(t0-t)/M ;
				double h = Math.pow( Math.sin(harg)/harg, M);
				double sarg = beta*(t0-t) ;
				double sin = Math.sin(sarg)/sarg;
				sum[0] += gAtBeta[i][0]*h*sin;
				sum[1] += gAtBeta[i][1]*h*sin;
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
