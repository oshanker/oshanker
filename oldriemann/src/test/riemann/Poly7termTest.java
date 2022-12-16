package riemann;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class Poly7termTest {
    static Poly7term poly7term012;
    static Poly7 poly7012;
    
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        poly7term012 = new Poly7term(0, 1, 2, 1, 0);
        poly7012 = new Poly7(0, 1, 2, 2, -1, 2);
    }
    
    @AfterClass
    public static void tearDownAfterClass() throws Exception {
    }
    
    @Test
    public void testDer() {
        double position = -0.1;
        double polyDer = poly7012.der(position);
        double f = poly7012.eval(position);
        Assert.assertEquals(2.63, polyDer, 1.0E-8);
        double termDer = poly7term012.der(position);
        double expected = f*f + 2*position*f*polyDer;
        Assert.assertEquals(expected, termDer, 1.0E-8);
    
        position = 0.5321546831790966;
        polyDer = poly7012.der(position);
        f = poly7012.eval(position);
        termDer = poly7term012.der(position);
        expected = f*f + 2*position*f*polyDer;
        Assert.assertEquals(expected, termDer, 1.0E-8);
        Assert.assertEquals(0, termDer, 1.0E-8);
    }
    
    @Test
    public void testEval() {
        double position = 1.05;
        double valplus = poly7term012.eval(position);
        double fplus = poly7012.eval(position);
        Assert.assertEquals(position*fplus*fplus, valplus, 1.0E-8);
    
        position = 0.5321546831790966;
        valplus = poly7term012.eval(position);
        fplus = poly7012.eval(position);
        Assert.assertEquals(position*fplus*fplus, valplus, 1.0E-8);
        double positionminus = 0.95;
        double valminus = poly7term012.eval(positionminus);
        double fminus = poly7012.eval(positionminus);
        Assert.assertEquals(positionminus*fminus*fminus, valminus, 1.0E-8);
    }
    
    @Test
    public void testPositionMax() {
        double positionMax05 = poly7term012.positionMax(0.5, 0, 1);
        Assert.assertEquals(0.5321546831790966, positionMax05, 1.0E-8);
//        double positionMax15 = poly7term012.positionMax(1.5, 1, 2);
//        Assert.assertEquals(1.577350269950018, positionMax15, 1.0E-8);
//        double max15 = poly7term012.eval(positionMax15);
//        Assert.assertEquals(-0.3849001794597505, max15, 1.0E-8);
        double max05 = poly7term012.eval(positionMax05);
        Assert.assertEquals(0.07106877367174713, max05, 1.0E-8);
//        System.out.println(max15 + " " + max05);
    }
}