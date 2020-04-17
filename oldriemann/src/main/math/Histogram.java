package math;


import java.text.NumberFormat;
import java.util.Random;

import riemann.NormalizedSpline;

public class Histogram {
    static final NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(3);
        nf.setMaximumFractionDigits(3);
        nf.setGroupingUsed(false);
    }
    
	final double min;
	final double max;
	final int binCount;
	final int[] hist;
	public final double delta;
	
	double sumSquares = 0;
	double sumY = 0;
	int sampleSize;
	
	public Histogram(double min, double max, int binCount) {
		this.min = min;
		this.max = max;
		this.binCount = binCount;
		hist = new int[binCount+2];
		delta = (max-min)/binCount;
	}

	public void addPoint(double y) {
		sampleSize++;
		hist[index(y)]++;
		sumSquares += y*y;
		sumY += y;
	}
	
	int index(double y) {
		if(y<min) {return 0;}
		if(y>=max) {return (hist.length-1);}
		int index = (int) ((y-min)/delta + 1);
		return index;
	}
	
	double yForIndex(int i) {
		if(i<=0) {return min-delta;}
		if(i>=hist.length-1) {return max;}
		double y = (i-1)*delta + min;
		return y;
	}
	
	
	/**
	 * x[i] = yForIndex(i+1)+delta/2;
	 * don't use end slots, since they are for overflow.
	 */
	public NormalizedSpline pdfSpline() {
		final double norm = sampleSize*delta;
		final double[] x = new double[binCount];
		final double[] y = new double[binCount];
		for(int i = 0; i < binCount; i++) {
			x[i] = yForIndex(i+1)+delta/2;
			y[i] = hist[i+1]/norm;
		}
        NormalizedSpline pdfSpline = new NormalizedSpline(x[0], x[binCount-1], y);	
		return pdfSpline;
	}
	
	/**
	 * x[i] = yForIndex(i)+delta/2;
	 */
	public double[] pdf() {
		final double norm = sampleSize*delta;
		
		final int length = hist.length;
		final double[] x = new double[length];
		final double[] y = new double[length];
		for(int i = 0; i < length; i++) {
			x[i] = yForIndex(i)+delta/2;
			y[i] = hist[i]/norm;
		}
		return y;
	}
	
	public double stdDev() {
		double mean = sumY/sampleSize;
		return Math.sqrt(sumSquares/sampleSize-mean*mean);
	}

	public double mean() {
		double mean = sumY/sampleSize;
		return mean;
	}

	public static void main(String[] args) {
		final double sigma = 1;
		final double mean = 0;
		final double min = mean - 4.5*sigma +0.25;
		final double max = mean + 4.5*sigma+0.25;
		final int binCount = 18;
		final Histogram hist = normalHist(sigma, mean, min, max, binCount);
		
		System.out.println(" mean " + nf.format(hist.mean()) + ", sigma "
				+ nf.format(hist.stdDev()));
		
		final double norm = hist.sampleSize*hist.delta;
		final double predNorm = 1.0/Math.sqrt(2*Math.PI*sigma);
		double sum = 0;
		double prevIncrement = 0.0d;
		double increment = 0.0d;
		
		final int length = hist.hist.length;
		for(int i = 0; i < length; i++) {
			increment = hist.hist[i];
			if(i==0) {
				sum = increment;
			} else if(i==length-1) {
				sum += prevIncrement/2.0 + increment;			
			} else if(i==1) {
				sum += increment/2.0;			
			} else {
				sum += prevIncrement/2.0 + increment/2.0;
			}
			prevIncrement = increment;
			final double x = hist.yForIndex(i)+hist.delta/2;
			final double pred = Math.exp(-(x-mean)*(x-mean)/(2*sigma))*predNorm;
			System.out.println(nf.format(x) + ", " 
			  + nf.format((double)hist.hist[i]/norm) 
			  + ", " + nf.format(pred)
			  + ", " + nf.format(((double)sum)/hist.sampleSize));
		}		

	}

	public static Histogram normalHist(final double sigma, final double mean, 
			final double min, final double max,
			final int binCount) {
		final int sampleCount = 1000000;
		
		return normalHist(sigma, mean, min, max, binCount, sampleCount);
	}

	public static Histogram normalHist(final double sigma, final double mean, final double min, final double max,
			final int binCount, final int sampleCount) {
		final Histogram hist = new Histogram(min, max, binCount);
		
		final Random generator = new Random(1781);
		for(int i = 0; i < sampleCount; i++) {
			final double y = sigma*generator.nextGaussian() + mean;
			hist.addPoint(y);
		}
		return hist;
	}

}
