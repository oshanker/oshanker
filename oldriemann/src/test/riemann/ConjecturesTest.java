package riemann;

import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.Random;

import org.junit.Test;

public class ConjecturesTest {

	@Test
	public void testCalculateDistribution() {
		Random x = new Random(5);
		int size = 1000002;
		char[] series = new char[size];
		int[][] doubletcounts = new int[2][4];
		int[][] tripletcounts = new int[2][8];
		char old1 = '1';
		char old = '1';
		for (int i = 0; i < series.length; i++) {
			float y = x.nextFloat();
			char current = '1';
			if(y<0.2){
				current = i%2==0?'-':'+';
			} else {
				current = i%2==0?'+':'-';
			}
			series[i] = current;
			if(i==0){
				old = current;
				continue;
			}
			int idx = 0;
			switch (old) {
			case '-':
				idx = 0;
				break;
			case '+':
				idx = 2;
				break;

			default:
				throw new IllegalStateException();
			}
			switch (current) {
			case '-':
				idx += 0;
				break;
			case '+':
				idx += 1;
				break;

			default:
				throw new IllegalStateException();
			}
			doubletcounts[(i-1)%2][idx]++;
			if(i==1){
				old1 = old;
				old = current;
				continue;
			}
			switch (old1) {
			case '-':
				idx += 0;
				break;
			case '+':
				idx += 4;
				break;

			default:
				throw new IllegalStateException();
			}
			tripletcounts[(i)%2][idx]++;
			old1 = old;
			old = current;
		}
		//System.out.println(Arrays.toString(series));
		System.out.println(Arrays.toString(Conjectures.descriptions[1]));
		for (int i = 0; i < doubletcounts.length; i++) {
			System.out.println(Arrays.toString(doubletcounts[(i+1)%2]));
		}
		System.out.println(Arrays.toString(Conjectures.descriptions[2]));
		for (int i = 0; i < tripletcounts.length; i++) {
			System.out.println(Arrays.toString(tripletcounts[(i+1)%2]));
		}
	}

}
