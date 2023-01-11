package riemann;

import java.text.NumberFormat;

public class SharpTransition {
    static NumberFormat nf = NumberFormat.getInstance();
    static {
        nf.setMinimumFractionDigits(7);
        nf.setMaximumFractionDigits(7);
        nf.setGroupingUsed(false);
    }

    public static void main(String[] args) throws Exception {
        double[] phiValues = {-0.2, -0.1, 0.0, 0.1, 0.2,};
        double[][] ratio = {
            {2 , 2.268576, 1.504897, 0.999090, 0.663992, 0.441400 , },
            {3 , 3.624520, 1.896329, 0.998864, 0.526210, 0.274796 , },
            {4 , 5.588850, 2.367189, 1.001220, 0.426549, 0.178314 , },
            {5 , 8.849518, 2.923269, 1.011907, 0.343228, 0.115100 , },
            {6 , 14.373004, 3.728921, 1.008597, 0.266224, 0.070801 , },
            {7 , 23.961623, 4.721631, 0.974761, 0.205256, 0.041497 , },
            {8 , 51.790514, 6.714706, 0.983726, 0.149338, 0.018210 , },
            {9 , 116.632353, 10.129730, 0.975052, 0.103749, 0.008631 , },
            {10 , 361.615385, 13.306122, 0.949264, 0.057506, 0.004515 , },
            {11 , 1397.000000, 34.533333, 1.027888, 0.041000, 0.001084 , },
        };
        for (int i = 0; i < ratio.length; i++) {
            int l = (int) ratio[i][0];
            System.out.print(l + " ");
            for (int j = 1; j < ratio[0].length; j++) {
                double f = ratio[i][j];
                double phi = phiValues[j-1];
                //double c = -Math.log(f)/(l*Math.tan(phi*Math.PI));
                //double c = -Math.log(f)/(l*Math.tan(phi));
                //y = 0.0383x + 1.9678
                double c = -Math.log(f)/(l*Math.tan(phi));
                System.out.print(nf.format(c) + " ");
            }
            System.out.println();
            
        }
    
    }
}
