package riemann;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class Poly7termTest {
    static Poly7term poly7term012;
    static Poly7 poly7012;
    static Poly7 poly7A012;
    static double A = 0.5;
    
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        poly7term012 = new Poly7term(0, 1, 2, 1, 5, 5);
        poly7012 = new Poly7(0, 1, 2, 2, -1, 2);
        poly7A012 = new Poly7(0, 1, 2, 2, -1, 2);
        poly7A012.setTerm(A, A*5, 5);
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
    
        position = 1.610702458700179;
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
    
        position = 1.610702458700179;
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
        double positionMax15 = poly7term012.positionMax(1.5, 1, 2);
        Assert.assertEquals(1.610702458700179, positionMax15, 1.0E-8);
        double max15 = poly7term012.eval(positionMax15);
        Assert.assertEquals(0.23619387078033371, max15, 1.0E-8);
        double max05 = poly7term012.eval(positionMax05);
        Assert.assertEquals(0.07106877367174713, max05, 1.0E-8);
    }
    
    @Test
    public void testPositionMaxA() {
        double positionMax05 = poly7A012.positionMax(0.5, 0, 1);
        Assert.assertEquals(0.44117451458624335, positionMax05, 1.0E-8);
        double positionMax15 = poly7A012.positionMax(1.5, 1, 2);
        Assert.assertEquals(1.5252718729923782, positionMax15, 1.0E-8);
        double max05 = poly7A012.eval(positionMax05);
        double expected0 = 0.41689197105413617;
        Assert.assertEquals(expected0, max05, 1.0E-8);
        double max15 = poly7A012.eval(positionMax15);
        double expected1= -0.2700198241673988;
        Assert.assertEquals(expected1, max15, 1.0E-8);
    }
}