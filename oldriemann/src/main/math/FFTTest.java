package math;
import org.jtransforms.fft.DoubleFFT_1D;

/**
 * Exercise some FFT code.
 * @author oshanker
 *
 */
public class FFTTest {
    public static void main(String[] args) {
        double[] input = new double[]{
                0.0176,
                -0.0620,
                0.2467,
                0.4599,
                -0.0582,
                0.4694,
                0.0001,
                -0.2873};
        DoubleFFT_1D fftDo = new DoubleFFT_1D(input.length);
        double[] fft = new double[input.length ];
        System.arraycopy(input, 0, fft, 0, input.length);
        fftDo.realForward(fft);

        for(double d: fft) {
            System.out.println(d);
        }
		System.out.println("*** get back ***");
		fftDo.realInverse(fft,true);
		for (double d: fft)
		System.out.println(d);
    }
}