package riemann;

import static org.junit.Assert.assertEquals;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class Poly7Test  {
    static Poly7 poly7012;
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        poly7012 = new Poly7(0, 1, 2, 2, -1, 2);
    }
    
    @AfterClass
    public static void tearDownAfterClass() throws Exception {
    }
    
    @Test
    public void testSecondDer(){
        double t = 1.5;
        double second = 6*t- 2*(poly7012.a + poly7012.b + poly7012.c);
        Assert.assertEquals(second , poly7012.secondDer(t), 1.0E-8);
    }
    
    @Test
    public void testDer() {
        double x2coeff = poly7012.t0 + poly7012.t1 + poly7012.t2;
        Assert.assertEquals(0, x2coeff, 1.0E-8);
        double x1coeff = poly7012.t0* (poly7012.b+ poly7012.c) +
            poly7012.t1* (poly7012.a+ poly7012.c) + poly7012.t2*(poly7012.a+ poly7012.b);
        Assert.assertEquals(0, x1coeff, 1.0E-8);
        double x0coeff = poly7012.t0* poly7012.b* poly7012.c +
            poly7012.t1* poly7012.a* poly7012.c + poly7012.t2* poly7012.a* poly7012.b;
        Assert.assertEquals(1, x0coeff, 1.0E-8);
        
        Assert.assertEquals(2.63, poly7012.der(-0.1), 1.0E-8);
    }
    
    @Test
    public void testEval() {
        double position = 1.05;
        double valplus = poly7012.eval(position);
        Assert.assertEquals(-0.0498750, valplus, 1.0E-8);
        double positionminus = 0.95;
        double valminus = poly7012.eval(positionminus);
        Assert.assertEquals(0.0498750, valminus, 1.0E-8);
    }
    
    @Test
    public void testFit() {
		Poly7 poly7 = new Poly7(
    			0.8749480149,
    			1.15890691298068,
    			1.3675025848634, 17.619276585379914,
    			-20.007604626096071598,
    			19.343950349024609636
        );
        poly7.setExtrema(    
        		 1.9266754104451154,
        	     -1.232146174810101691,      0

        );
        
        double deviation = poly7.setTermValues();
        assertEquals("deviation large", 0.0 , deviation, 1.0E-8);
        
        double t1 = 244.02115917156453;
        double incr = 0.24359904690399015;
        double[] zeta = {
        		1.926498073039712
        		, -1.23155062255226};

        for (int i = 0; i < 2; i++) {
			double z = poly7.eval(t1-243);
            System.out.println(t1 + " zetaFromRiemann " + 
              zeta[i] + " poly7 " + z);
            assertEquals("not equal", zeta[i] , z, 1.0E-4);
            t1+=incr;
         }
        
       
        System.out.println("============= " );
    	
    }
    
    @Test
    public void testPositionMax() {
        double positionMax05 = poly7012.positionMax(0.5, 0, 1);
        Assert.assertEquals(0.42264972993696365, positionMax05, 5.0E-7);
        double positionMax15 = poly7012.positionMax(1.5, 1, 2);
        Assert.assertEquals(1.577350269950018, positionMax15, 2.0E-8);
        Assert.assertEquals(2.0, (positionMax05+positionMax15), 1.0E-8);
        double max15 = poly7012.eval(positionMax15);
        Assert.assertEquals(-0.3849001794597505, max15, 1.0E-8);
        double max05 = poly7012.eval(positionMax05);
        Assert.assertEquals(0.3849001794597505, max05, 1.0E-8);
    }
}