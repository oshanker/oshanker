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
import java.util.List;

/**
 * This is an example on how to extract text line by line from pdf document
 */
public class GetLinesFromPDF extends PDFTextStripper {
	
	List<String> lines = new ArrayList<String>();
	int count = 0;

	public GetLinesFromPDF() throws IOException {
	}

	/**
	 * @throws IOException If there is an error parsing the document.
	 */
	public static void main( String[] args ) throws IOException	{
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
            System.out.println("getNumberOfPages " +  document.getNumberOfPages() );
            for (int i = 2; i < 4; i++) {
                stripper.setStartPage( i );
                stripper.setEndPage( i);

                Writer dummy = new OutputStreamWriter(new ByteArrayOutputStream());
                stripper.writeText(document, dummy);
                
                // print lines
                for(String line:stripper.lines){
                    System.out.println(line);               
                }
                stripper.lines.clear();
                
                System.out.println("*********************************");
            }
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
		lines.add(str);
		// you may process the line here itself, as and when it is obtained
	}
}