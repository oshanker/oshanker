package math;

import java.text.NumberFormat;
import java.util.Random;

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
	private final double delta;
	
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
	
	public double stdDev() {
		double mean = sumY/sampleSize;
		return Math.sqrt(sumSquares/sampleSize-mean*mean);
	}

	public double mean() {
		double mean = sumY/sampleSize;
		return mean;
	}

	public static void main(String[] args) {
		double sigma = 1.7;
		double mean = 2.3;
		double min = mean - 3*sigma;
		double max = mean + 3*sigma;
		int binCount = 6;
		Histogram hist = new Histogram(min, max, binCount);
		for(double y = min-sigma; y<= max; y += sigma) {
			System.out.println(nf.format(y) + "," + hist.index(y));
		}
		
		Random generator = new Random(1781);
		double sumSquares = 0;
		int sampleCount = 1000000;
		for(int i = 0; i < sampleCount; i++) {
			double y = sigma*generator.nextGaussian() + mean;
			sumSquares += y*y;
			hist.addPoint(y);
		}
		
		System.out.println(" index of mean " + hist.index(mean) + ", sigma "
				+ nf.format(Math.sqrt(sumSquares/hist.sampleSize-mean*mean)));
		System.out.println(" mean " + nf.format(hist.mean()) + ", sigma "
				+ nf.format(hist.stdDev()));
		
		int sum = 0;
		for(int i = 0; i < hist.hist.length; i++) {
			sum += hist.hist[i];
			System.out.println(nf.format(hist.yForIndex(i)) + ", " + hist.hist[i]+ ", " + sum);
		}		

	}

}
