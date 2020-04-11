package math;

import static org.junit.Assert.*;

import org.junit.Test;

import riemann.NormalizedSpline;

public class HistogramTest {

	@Test
	public void testStdDev() {
		final double sigma = 1;
		final double mean = 0;
		final double min = mean - 4.5*sigma +0.125;
		final double max = mean + 4.5*sigma+0.125;
		final int binCount = 36;
		final int sampleCount = 2000000;
		final Histogram hist = Histogram.normalHist(sigma, mean, min, max, binCount, sampleCount);
		
		
		final double norm = hist.sampleSize*hist.delta;
		final double predNorm = 1.0/Math.sqrt(2*Math.PI*sigma);
		
		final int length = hist.hist.length;
		final double[] x = new double[length];
		double[] y = new double[length];
		double error = 0;
		for(int i = 0; i < length; i++) {
			x[i] = hist.yForIndex(i)+hist.delta/2;
			final double pred = Math.exp(-(x[i]-mean)*(x[i]-mean)/(2*sigma))*predNorm;
			y[i] = hist.hist[i]/norm;
			assertEquals(pred, y[i], 0.0011);
			double currentError = Math.abs(pred - y[i]);
			if(error<currentError) {error = currentError;}
			if(Math.abs(x[i]) <= 0.25) {
				System.out.println(Histogram.nf.format(x[i]) + ", " 
			  + Histogram.nf.format(y[i]));
			}
		}	
		System.out.println(error + ", hist.delta " + hist.delta);
		
        NormalizedSpline normalizedSpline = new NormalizedSpline(x[0], x[length-1], y);	
		for(int i = 1; i < length-2; i++) {
			double xi = hist.yForIndex(i)+ 0.75*hist.delta;
			final double pred = Math.exp(-(xi-mean)*(xi-mean)/(2*sigma))*predNorm;
			final double yi = normalizedSpline.eval(xi);
			assertEquals(pred, yi, 0.0012);
		}	
        
    }

}
