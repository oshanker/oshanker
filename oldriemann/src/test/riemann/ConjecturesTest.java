package riemann;

import static org.junit.Assert.*;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.Random;

import org.junit.Test;

import math.GSeries;

public class ConjecturesTest {

    @Test
    public void testX() {
        double t = 1.0E12 + 244.158906912980683962;
        increment(t);
        double incr = increment(1.0E15 + 192.309350419702134727 );
        BigDecimal tval = new BigDecimal(192.2043554309546, Gram.mc).add(
                BigDecimal.valueOf(1.0E15), Gram.mc);
        BigDecimal t2 = tval.divide(Gram.bdTWO);
        BigDecimal sqrtArg1 = Gram.sqrt(tval.divide(Gram.pi_2, Gram.mc), Gram.mc, 1.0E-25);
        double basesqrtArg1 = sqrtArg1.doubleValue();
        int k1 = (int)basesqrtArg1;
        Gram.initLogVals(k1);
        BigDecimal lnsqrtArg1BD = Gram.log(sqrtArg1, Gram.mc);
        BigDecimal thetaPi = tval.multiply(lnsqrtArg1BD, Gram.mc).subtract(t2, Gram.mc)
                .subtract(Gram.pi8, Gram.mc).divide(Gram.pi, Gram.mc);
        System.out.println(thetaPi);
        /*
baseGram 192.2043554309546 idx 5045354828590534
        double[] begin = {192.309350419702134727, 1921602.793342093316678265};
f  : [1.9524278734805527, 0.2869931460545315] theta 1.716715278515568 rotatedSum 1.1379834185443194E-4 zeta 7.381699365530212E-11
zetaFromRiemann 7.413854952046176E-11
sqrtArg1[i].doubleValue() 1.2615662610102013E7 correction -1.1379826803743828E-4
g  : [-1.23276504691923, 0.32595788950860927] argalphaBase 2.4139476198307377
f  : [0.7037596046854905, -1.063334672364068] theta 0.5846098679289753 rotatedSum 1.1521610619680267E-4 zeta -1.3689409652075171E-8
zetaFromRiemann 4.895887538723708E-9
sqrtArg1[i].doubleValue() 1.2615662622221947E7 correction -1.1522979560645475E-4
         */
    }

    private double increment(double t) {
        double incr = Math.PI/Math.log(Math.sqrt(t/(2*Math.PI)));
        System.out.println(incr);
        t += 1000000*incr;
        incr = Math.PI/Math.log(Math.sqrt(t/(2*Math.PI)));
        System.out.println(incr);
        return incr;
    }

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
