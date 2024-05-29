package riemann;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

public class Poly7termBTest {
    static Poly7term poly7term012;
    static Poly7 poly7012;
    static Poly7 poly7B012;
    static double B = 0.5;
    @BeforeClass
    public static void setUpBeforeClass() throws Exception {
        poly7term012 = new Poly7term(0, 1, 2, 0, 1, 0);
        poly7012 = new Poly7(0, 1, 2, 2, -1, 2);
        poly7B012 = new Poly7(0, 1, 2, 2, -1, 2);
        poly7B012.setTerm(0, B, 0);
    }
    
    @AfterClass
    public static void tearDownAfterClass() throws Exception {
    }
    
    @Test
    public void testDer() {
        double position = -0.1;
        double polyDer = poly7012.der(position);
        double f = poly7012.eval(position);
        double termDer = poly7term012.der(position);
        double expected =  2*f*polyDer;
        Assert.assertEquals(expected, termDer, 1.0E-8);
    
        position = 0.42264973081044926;
        polyDer = poly7012.der(position);
        f = poly7012.eval(position);
        termDer = poly7term012.der(position);
        expected = 2*f*polyDer;
        Assert.assertEquals(expected, termDer, 1.0E-8);
        Assert.assertEquals(0, termDer, 1.0E-8);
    
        position = 1.5773502691895684;
        polyDer = poly7012.der(position);
        f = poly7012.eval(position);
        termDer = poly7term012.der(position);
        expected =  2*f*polyDer;
        Assert.assertEquals(expected, termDer, 1.0E-8);
        Assert.assertEquals(0, termDer, 1.0E-8);
    }
    
    @Test
    public void testEval() {
        double position = 1.05;
        double valplus = poly7term012.eval(position);
        double fplus = poly7012.eval(position);
        Assert.assertEquals(fplus*fplus, valplus, 1.0E-8);
    
        position = 0.42264973081044926;
        valplus = poly7term012.eval(position);
        fplus = poly7012.eval(position);
        Assert.assertEquals(fplus*fplus, valplus, 1.0E-8);
    
        position = 1.5773502691895684;
        valplus = poly7term012.eval(position);
        fplus = poly7012.eval(position);
        Assert.assertEquals(fplus*fplus, valplus, 1.0E-8);
    
        double positionminus = 0.95;
        double valminus = poly7term012.eval(positionminus);
        double fminus = poly7012.eval(positionminus);
        Assert.assertEquals(fminus*fminus, valminus, 1.0E-8);
    }
    
    @Test
    public void testPositionMax() {
        double positionMax05 = poly7term012.positionMax(0.5, 0, 1);
        Assert.assertEquals(0.42264973081044926, positionMax05, 1.0E-8);
        double positionMax15 = poly7term012.positionMax(1.5, 1, 2);
        Assert.assertEquals(1.5773502691895684, positionMax15, 1.0E-8);
        double max15 = poly7term012.eval(positionMax15);
        Assert.assertEquals(0.1481481481481481, max15, 1.0E-8);
        double max05 = poly7term012.eval(positionMax05);
        Assert.assertEquals(0.1481481481481481, max05, 1.0E-8);
    }
    
    @Test
    public void testPositionMaxB() {
        double positionMax05 = poly7B012.positionMax(0.5, 0, 1);
        Assert.assertEquals(0.42264973081044926, positionMax05, 1.0E-7);
        double positionMax15 = poly7B012.positionMax(1.5, 1, 2);
        Assert.assertEquals(1.5773502691895684, positionMax15, 1.0E-7);
        double max05 = poly7B012.eval(positionMax05);
        double expected0 = 0.3849001794597505 + B*0.1481481481481481;
        Assert.assertEquals(expected0, max05, 1.0E-7);
        double max15 = poly7B012.eval(positionMax15);
        double expected1= -0.3849001794597505 + B*0.1481481481481481;
        Assert.assertEquals(expected1, max15, 1.0E-7);
    }
}