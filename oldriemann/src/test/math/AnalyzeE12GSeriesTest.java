package math;

import org.junit.Assert;
import org.junit.Test;
import riemann.CopyZeroInformation;
import riemann.Interpolate;
import riemann.Poly4;
import riemann.Rosser;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.LinkedList;

import static math.AnalyzeE12GSeries.testGetSavedGSeries1;
import static math.FixE12GSeries.TEST_VALES;

public class AnalyzeE12GSeriesTest  {
    @Test
    public void testFixGSeriesUsingPosmax() {
        GSeries gAtBeta = Interpolate.readGSeries();
        String zerosFile = Rosser.getParam("zerosFile");
        BufferedReader[] zeroIn = AnalyzeE12GSeries.zerosFileWithMaxPos(zerosFile);
        AnalyzeE12GSeries.fixGSeriesUsingPosmax(zeroIn, gAtBeta);
    }
    
    @Test
    public void testChangeToZetaAndDer() {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        for (int i = 0; i < TEST_VALES.length; i++) {
            zeroInfo.add(TEST_VALES[i]);
        }
        //AnalyzeE12GSeries analyzeE12GSeries = new AnalyzeE12GSeries();
        GSeries gAtBeta = Interpolate.readGSeries();
        double[] deviation1 = AnalyzeE12GSeries.testChangeToZetaAndDer(
            zeroInfo, 1999907, gAtBeta, false);
        Assert.assertTrue(deviation1[0] < 1.0E-9);
        System.out.println(" " + Arrays.toString(deviation1));
    }
    
    @Test
    public void testChangeToZetaAndDerBad() {
        double[][] badZeros = {
            {5195.318244135411, 7.516433369861296, 0.9165390902393957, 5195.478496162874},
            {5195.672422220875, -1.3049505772938506, -0.009900533082844082, 5195.687646467244},
            {5195.703116888218, 1.2633640076170625, 0.4813896177508429, 5195.853812920017},
        };
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        for (int i = 0; i < badZeros.length; i++) {
            zeroInfo.add(badZeros[i]);
        }
        GSeries gAtBeta = Interpolate.readGSeries();
        double[] deviation1 = AnalyzeE12GSeries.testChangeToZetaAndDer(
            zeroInfo, 40650, gAtBeta, true);
        System.out.println(" " + Arrays.toString(deviation1));
    }
    
    @Test
    public void testTestGetSavedGSeries1() {
        double firstZero = 248;
        final double stopValue = 243839.0;
        boolean ignoreMax = true;
        GSeries gAtBeta = Interpolate.readGSeries();
        String zerosFile = Rosser.getParam("zerosFile");
        BufferedReader[] zeroIn = AnalyzeE12GSeries.zerosFileWithMaxPos(zerosFile);
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        testGetSavedGSeries1(firstZero, zeroIn, gAtBeta,
            2000002, stopValue, ignoreMax, zeroInfo);
        System.out.println(Arrays.toString(zeroInfo.get(0)));
    }
    
    @Test
    public void savePositionMax() throws IOException {
        LinkedList<double[]> zeroInfo = new LinkedList<>();
        FileOutputStream fileOutputWriter = new FileOutputStream("out/gSeriesE12/positionMax.csv");
        PrintStream outputStream = new PrintStream(fileOutputWriter);
        for (int i = 0; i < 1000004; i++) {
            double[] nextValues = CopyZeroInformation.skipUntil(Interpolate.zeroIn, 0);
            if (nextValues == null) {
                break;
            }
            if(nextValues[1]<0){
                nextValues[2]=-nextValues[2];
            }
            if (i>0) {
                double[] oldZero = zeroInfo.getLast();
                Poly4 poly = new Poly4(oldZero[0], nextValues[0], oldZero[1], nextValues[1],
                    oldZero[2]);
                double positionMax = poly.getPositionMax();
                outputStream.println(positionMax);
            }
            zeroInfo.add(nextValues);
        }
        outputStream.close();
    }
}