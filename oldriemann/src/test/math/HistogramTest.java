package math;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import riemann.NormalizedSpline;

public class HistogramTest {
	final double sigma = 1;
	final double mean = 0;
	final double min = mean - 4.5*sigma +0.125;
	final double max = mean + 4.5*sigma+0.125;
	final int binCount = 36;
	final int sampleCount = 2000000;
	final Histogram hist = Histogram.normalHist(sigma, mean, min, max, binCount, sampleCount);
	final int length = hist.hist.length;
    final double epsilon = 0.0005;
	
	
	final double predNorm = 1.0/Math.sqrt(2*Math.PI*sigma);

	@Test @Ignore
	public void testStdDev() {
		
		final double[] x = new double[length];
		final double[] y = hist.pdf();
		double error = 0;
		for(int i = 0; i < length; i++) {
			x[i] = hist.yForIndex(i)+hist.delta/2;
			final double pred = Math.exp(-(x[i]-mean)*(x[i]-mean)/(2*sigma))*predNorm;
			assertEquals(pred, y[i], 0.0011);
			double currentError = Math.abs(pred - y[i]);
			if(error<currentError) {error = currentError;}
			if(Math.abs(x[i]) <= 0.25) {
				System.out.println(Histogram.nf.format(x[i]) + ", " 
			  + Histogram.nf.format(y[i]));
			}
		}	
		System.out.println(error + ", hist.delta " + hist.delta);
		
        NormalizedSpline pdfSpline = hist.pdfSpline();	
		for(int i = 2; i < length-3; i++) {
			double xi = hist.yForIndex(i)+ 0.75*hist.delta;
			final double pred = Math.exp(-(xi-mean)*(xi-mean)/(2*sigma))*predNorm;
			final double yi = pdfSpline.eval(xi);
			assertEquals(pred, yi, 0.0012);
		}	
        
    }


	@Test(expected = IllegalArgumentException.class)
	public void testCDFException() {
		hist.findQuartile(1.1, epsilon);
	}
	
	@Test
	public void testCDF() {
        NormalizedSpline cdfSpline = hist.cdfSpline();	
		final double xlow = hist.yForIndex(1)+ 0.5*hist.delta;
		final double xhigh = hist.yForIndex(hist.hist.length-2);
		double x = cdfSpline.findX(xlow, xhigh, 0.025, epsilon );
		System.out.println(x);
		assertEquals(-1.96, x, 0.02);
		x = hist.findQuartile(0.025, epsilon);
		assertEquals(-1.96, x, 0.02);
		System.out.println(x);
		
		double xmid = hist.findQuartile(0.5, epsilon);
		System.out.println(xmid);
		assertEquals(0, xmid, 0.025);
		
		x = cdfSpline.findX(xlow, 
				xhigh, 0.975, epsilon );
		System.out.println(x);
		assertEquals(1.96, x, 0.025);
		x = hist.findQuartile(0.975, xmid, epsilon);
		System.out.println(x);
		assertEquals(1.96, x, 0.025);
		
	
	}
}
