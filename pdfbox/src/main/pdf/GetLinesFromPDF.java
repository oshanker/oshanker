package pdf;

import org.apache.pdfbox.pdmodel.PDDocument;
import org.apache.pdfbox.text.PDFTextStripper;
import org.apache.pdfbox.text.TextPosition;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * This is an example on how to extract text line by line from pdf document
 */
public class GetLinesFromPDF extends PDFTextStripper {
    
    int count = 0;
    int[] hist = new int[10];
    static enum TYPE{k,t,Z,Other, kHeader, tHeader};
    TYPE currentType = TYPE.Other;

    public GetLinesFromPDF() throws IOException {
    }

    /**
     * @throws IOException If there is an error parsing the document.
     */
    public static void main( String[] args ) throws IOException {
        PDDocument document = null;
        String fileName = "data/peakZ15-16.pdf";
        File file = new File(fileName);
        if (!file.exists()) {
            System.out.println("Doesn't exist:" + file.getAbsolutePath());
            System.exit(-1);
        }

        try {
            document = PDDocument.load( new File(fileName) );
            GetLinesFromPDF stripper = new GetLinesFromPDF();
            stripper.setSortByPosition( true );
            int numberOfPages = document.getNumberOfPages();
            System.out.println("getNumberOfPages " +  numberOfPages );
            for (int i = 2; i <= numberOfPages; i++) {
                stripper.setStartPage( i );
                stripper.setEndPage( i);
                stripper.currentType = TYPE.Other;

                Writer dummy = new OutputStreamWriter(new ByteArrayOutputStream());
                stripper.writeText(document, dummy);
                
            }
            System.out.println(stripper.count + ", " + Arrays.toString(stripper.hist));
        }
        finally {
            if( document != null ) {
                document.close();
            }
        }
    }

    /**
     * Override the default functionality of PDFTextStripper.writeString()
     */
    @Override
    protected void writeString(String str, List<TextPosition> textPositions) throws IOException {
        Pattern p = Pattern.compile("[0-9].*");
        Matcher m = p.matcher(str);
        switch (currentType) {
        case k:
            if(m.matches()){
               currentType = TYPE.t;
            }
            return;
        case t:
            currentType = TYPE.Z;
            return;
        case Z:
            currentType = TYPE.k;
            count++;
            double Z = Math.abs(Double.parseDouble(str));
            int idx = (int) ((Z-500)/50);
            if(idx >= hist.length){idx = hist.length - 1; }
            hist[idx]++;
            break;
        case kHeader:
            if(str.equals("t")){
               currentType = TYPE.tHeader;
            } else {
                currentType = TYPE.Other;
            }
            break;
        case tHeader:
            if(str.startsWith("Z(t)")){
               currentType = TYPE.k;
            } else {
                currentType = TYPE.Other;
            }
            break;

        default:
            if(str.equals("k")){
                currentType = TYPE.kHeader;
            }
            return;
        }
        // you may process the line here itself, as and when it is obtained
    }
}