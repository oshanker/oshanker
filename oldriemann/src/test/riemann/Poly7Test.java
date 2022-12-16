package riemann;

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
    public void testPositionMax() {
        double positionMax05 = poly7012.positionMax(0.5, 0, 1);
        Assert.assertEquals(0.42264972993696365, positionMax05, 1.0E-8);
        double positionMax15 = poly7012.positionMax(1.5, 1, 2);
        Assert.assertEquals(1.577350269950018, positionMax15, 1.0E-8);
        Assert.assertEquals(2.0, (positionMax05+positionMax15), 1.0E-8);
        double max15 = poly7012.eval(positionMax15);
        Assert.assertEquals(-0.3849001794597505, max15, 1.0E-8);
        double max05 = poly7012.eval(positionMax05);
        Assert.assertEquals(0.3849001794597505, max05, 1.0E-8);
        System.out.println(max15 + " " + max05);
    }
}