# read example

https://github.com/oshanker/oshanker/commit/660c141ce2768e585b3926ad29f67a8001ae4caa#diff-07e2157e60a134ec1753b14713faba1c556a860c27ad8d846647292e4b5dec6c 

oldriemann/src/test/math/GSeriesTest.java 
~~~java
   private void testReadE12(int sampleIndex, PrintWriter out) throws Exception{
        double t0 = gramE12[sampleIndex][0];
        int index = (int) Math.floor(t0);
        BigDecimal offset = BigDecimal.valueOf(1.0E12);
        int k0 = 1, k1=398942;
        File file = new File("data/gSeriesE12/" + Integer.toString(index) +".dat");
        InputStream is = new FileInputStream(file);
        // create buffered input stream.
        BufferedInputStream bis = new BufferedInputStream(is);

        // create data input stream to read data in form of primitives.
        DataInputStream in = new DataInputStream(bis);
        final int initialPadding = 40;
        int R = 30000+2*initialPadding;
        double begin = in.readDouble();
        double incr = in.readDouble();
        double[][] gBeta = new double[R][2];
        for (int i = 0; i < gBeta.length; i++) {
            gBeta[i][0] = in.readDouble();
            gBeta[i][1] = in.readDouble();
        }
        GSeries gAtBeta = new GSeries(k0, k1, offset,  begin,  incr, gBeta);
        in.close();
~~~
