package riemann;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Deque;
import java.util.LinkedList;

import org.junit.Ignore;
import org.junit.Test;

public class RosserTest {

    @Test @Ignore
    public void testThreeZeros() throws IOException {
        Rosser.readConfig("data/RosserConfig.txt");
        BufferedReader[] zeroIn = Rosser.getZerosFile();
        String input = zeroIn[0].readLine();
        double zero = getZeroFromFile(input);
        System.out.printf("first zero %f ", zero);
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        System.out.printf("gram %f \n", (baseLimit -zero)/gramIncr);
        int gramOffset = (int) Math.floor((baseLimit - zero) / gramIncr);
        double beginGram = baseLimit - gramOffset*gramIncr;
        System.out.printf("beginGram %f, check  %f\n", beginGram, beginGram + gramOffset*gramIncr);
        while (zero < beginGram) {
            input = zeroIn[0].readLine();
            zero =getZeroFromFile(input);
            System.out.printf("next zero %f \n", zero);
        }
        PrintStream out = getOutputPrintStream();
        Deque<Integer> buffer = new LinkedList<>();

        int count = 1; //244.158907
        int indexOfThree = 1000;
        double nextGram = beginGram + gramIncr;
            for (int i = 0; i < 10000000; i++) {
                while (zero < nextGram) {
                    String input1 = zeroIn[0].readLine();

                    zero = getZeroFromFile(input1);
                    if(zero==-1){
                        break;
                    }
                    if (zero >= nextGram) {
                        break;
                    }
                    count++;
                }
                indexOfThree = addAndResize(buffer, count, indexOfThree);
                //((LinkedList<Integer>) buffer).get(19);
                if(buffer.size() == 29 && indexOfThree <= 19) {
                    dumpBufferAndClear(out, (LinkedList<Integer>) buffer);
                    indexOfThree = 1000;
                }
                // handle empty
                beginGram = nextGram;
                nextGram = beginGram + gramIncr;
                count = 0;
                while (zero >= nextGram) {
                    //while
                    ++i;
                    indexOfThree = addAndResize(buffer, count, indexOfThree);
                    if(buffer.size() == 29 && indexOfThree <= 19) {
                        dumpBufferAndClear(out, buffer);
                        indexOfThree = 1000;
                    }
                    beginGram = nextGram;
                    nextGram = beginGram + gramIncr;
                }
                count = 1;
            }

        out.close();
    }

    private double getZeroFromFile(String input) {
        String[] parsed = input.trim().split("\\s+");
        double zero = Double.parseDouble(parsed[0]);
        if (zero < 0) {
            return zero;
        }
        if (parsed.length > 1) {
            zero += Double.parseDouble(parsed[1]);
        }
        return zero;
    }

    private int addAndResize(Deque<Integer> buffer, int count, int indexOfThree) {
        buffer.addLast(count);
        //if buffer size > 29, remove first
        if (buffer.size() > 29) {
            buffer.removeFirst();
            if(indexOfThree == 0){
                indexOfThree = 1000;
            }
            if(indexOfThree<1000){
                indexOfThree--;
            }
        }
        if (count >= 3 && (indexOfThree==1000)) {
            indexOfThree = buffer.size()-1;
        }
        return indexOfThree;
    }

    private void dumpBufferAndClear(PrintStream out, Deque<Integer> buffer) {
        dumpBuffer(out, (LinkedList<Integer>) buffer);
        while(buffer.size() > 10) {
            buffer.removeFirst();
        }
    }

    private void dumpBuffer(PrintStream out, LinkedList<Integer> buffer) {
        int i = 0;
        out.print(buffer.get(i) );
        for (i = 1; i < buffer.size(); i++) {
            out.print("," + buffer.get(i) );
        }
        out.println();
    }

    @Test @Ignore
    public void testGetZerosFile() throws IOException {
        Rosser.readConfig("data/RosserConfig.txt");
        BufferedReader[] zeroIn = Rosser.getZerosFile();
        String input = zeroIn[0].readLine();
        double zero =getZeroFromFile(input);
        System.out.printf("first zero %f ", zero);
        double baseLimit = Rosser.getParamDouble("baseLimit");
        double gramIncr = Rosser.getParamDouble("gramIncr");
        System.out.printf("gram %f \n", (baseLimit -zero)/gramIncr);
        int gramOffset = (int) Math.floor((baseLimit - zero) / gramIncr);
        double beginGram = baseLimit - gramOffset*gramIncr;
        System.out.printf("beginGram %f, check  %f\n", beginGram, beginGram + gramOffset*gramIncr);
        while (zero < beginGram) {
            input = zeroIn[0].readLine();
            zero = getZeroFromFile(input);
            System.out.printf("next zero %f \n", zero);
        }
        PrintStream out = getOutputPrintStream();

        int count = 1; //244.158907
        int cumulative = 0;
        double nextGram = beginGram + gramIncr;
        for (int i = 0; i < 25615; i++) {
            while (zero < nextGram) {
                input = zeroIn[0].readLine();
                zero = getZeroFromFile(input);
                if (zero >= nextGram) {
                    break;
                }
                count++;
            }
            if (i==0) {
                out.print(count );
            } else {
                out.print("," + count);
            }
            // handle empty
            beginGram = nextGram;
            nextGram = beginGram + gramIncr;
            cumulative += count;
            count = 0;
            while (zero >= nextGram) {
                //while
                ++i;
                out.print("," + count);
                beginGram = nextGram;
                nextGram = beginGram + gramIncr;
            }
            count = 1;
        }

        System.out.printf("cumulative %d", cumulative);
        out.close();

    }

    private PrintStream getOutputPrintStream() {
        PrintStream out = null;
        File file = new File("../python/out/intervalsTestE28Threes.csv");
        if (!file.exists()) {
            try {
                file.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        try {
            out = new PrintStream(file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return out;
    }

    @Test
    public void testTypeIIRatios() throws Exception{
        int displacementCount = 5;
        double[][] typeIIratios = new double[10][displacementCount ];
        String ratiosFile = "data/typeIIratios.txt";
        BufferedReader ratiosIn = new BufferedReader(new FileReader(ratiosFile));
        for (int j = 0; j < typeIIratios.length; j++) {
            String val = ratiosIn.readLine();
            val = val.replaceAll("\\\\", "");
            String[] valsIn = val.split("&");
            for (int displacement = 0; displacement < displacementCount; displacement++) {
                typeIIratios[j][displacement] = Double.parseDouble(valsIn[displacement + 1]);
            }
            for (int displacement = 0; displacement < displacementCount; displacement++) {
                double val1 = Math.log(typeIIratios[j][displacement]);
                double phi = -(displacement-displacementCount/2)*Math.PI/10;
                double val2 = val1;
                if(displacement != displacementCount/2){
                    val2 = Math.tan(phi)*(1+Math.cos(2*phi)/4.5)*(j+2-0.5)/1.5;
                }
                System.out.print((displacement==0?(j+2+" & "): " & ") 
                        + Conjectures.nf.format(val1/val2 ) );
//                System.out.print((displacement==0?(j+2+" & "): " & ") 
//                        + Conjectures.nf.format(
//                                typeIIratios[j][displacement]/Math.exp(val2)) );
            }
            System.out.println(" \\\\"); 
        }
        ratiosIn.close();
    }

}
