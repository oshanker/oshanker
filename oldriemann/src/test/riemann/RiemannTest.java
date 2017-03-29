package riemann;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class RiemannTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testRiemannDoubleLong() {
		double[] args = {1.87383225, 2.13054990,        
		        2.38726754,   2.64398519        };
		double[] vals = {4.023671,42.85232072, -4.52434707, 
				-2.84388877};
		long init = 0;
		for (int i = 0; i < args.length; i++) {
			double zeta = Riemann.riemann(args[i], 267653395647L);
			System.out.println(args[i] + ":" +zeta + ":" + Math.abs(zeta-vals[i]));
			assertTrue(i + ":" + Math.abs(zeta-vals[i]), Math.abs(zeta-vals[i])<0.0000005);
			if(i==0){init = System.currentTimeMillis();}
		}
		long end = System.currentTimeMillis();
		System.out.println("N " + Riemann.oldN + " calc for " + (args.length-1) + " " + (end - init) + "ms");
	}

}
