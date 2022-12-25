package math;

public class LinearEquation
{
    double [][] coefficients;
    double [][] values;
    int rowIndex[];

    public static void main(String args[]) {
        LinearEquation linearEquation = new LinearEquation(7);
        linearEquation.runInvert();
    }

    public void runInvert() {
        int n = values.length;
        //Matrix representation
        for(int i=0; i<n; i++)
        {
            for(int j=0; j<n; j++)
            {
                System.out.print(" " + coefficients[i][j]);
            }
            System.out.print("  =  " + values[i][0]);
            System.out.println();
        }

        //inverse of matrix mat[][]
        double inverted_mat[][] = invert();
        System.out.println("The inverse is: ");
        printMatrix(inverted_mat);
        //Multiplication of mat inverse and constants
        //multiplyByInverse(inverted_mat);

    }
    
    public static void printMatrix(double[][] matrix) {
        int length = matrix.length;
        for (int i = 0; i < length; ++i)
        {
            for (int j = 0; j < matrix[i].length; ++j)
            {
                System.out.print(matrix[i][j] + "  ");
            }
            System.out.println();
        }
    }
    
    public static double[][] transpose(double[][] matrix) {
        int length = Math.min(matrix.length, matrix[0].length);
        double[][] transpose = new double[matrix[0].length][matrix.length];
        for (int row = 0; row < length; ++row)
        {
            transpose[row][row] = matrix[row][row];
            for (int j = row+1; j < matrix[row].length; ++j)
            {
                transpose[row][j] = matrix[j][row];
                transpose[j][row] = matrix[row][j];
            }
        }
        return transpose;
    }

    private void multiplyByInverse(double[][] inverted_mat) {
        int n = values.length;
        double result[][] = new double[n][1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < 1; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    result[i][j] = result[i][j] + inverted_mat[i][k] * values[k][j];
                }
            }
        }
        System.out.println("The solution is:");
        for(int i=0; i<n; i++)
        {
            System.out.println(result[i][0] + " ");
        }
    }

    public LinearEquation(double[][] coefficients, double[][] values) {
        this.coefficients = coefficients;
        this.values = values;
    }

    public LinearEquation(double[][] coeff) {
        this.coefficients = new double[coeff.length][coeff[0].length];
        for(int row = 0; row < coefficients.length; row++)
        {
            for(int col = 0; col < coefficients[0].length; col ++)
            {
                this.coefficients[row][col ] = coeff[row][col ];
            }
        }
        gaussian(this.coefficients);
    }

    public LinearEquation(int n) {
        coefficients = new double[n][n];
        values = new double[n][1];
        for(int row = 0; row < n; row++)
        {
            for(int col = 0; col <= row; col ++)
            {
                coefficients[row][col ] = 1;
            }
            values[row][0] = row + 100001.5;
        }
        coefficients[0][n-2 ] += 800000;
        values[0][0] += 800000 - 100000;

        coefficients[1][n-2] += 700000;
        values[1][0] += 700000;

        coefficients[n-1][n-2] += 700000;
        values[n-1][0] += 700000;

        coefficients[n-1][n-1 ] += 20;
        values[n-1][0] += 41;
    }

    public double[][] invert()
    {

        // Transform the matrix into an upper triangle
        gaussian(coefficients);

        //printMatrix(coefficients);

        //System.out.println(Arrays.toString(index));

        //double[][] transformFromIdentity = populateTransformFromIdentity(values, n);
        //double[][] transformFromIdentity = populateValueHolder(values, n);

        double[][] transformFromIdentity = values;


        double[][] inverse = solveDoubleArray( transformFromIdentity);
        return inverse;
    }

    public static double[] multiply(double[][] matrix, double[] vector) {
        double[] ret = new double[matrix.length];
        for (int row = 0; row < matrix.length; row++) {
            ret[row] = 0;
            for (int col = 0; col < vector.length; col++) {
                ret[row] += matrix[row][col]*vector[col];
            }
        }
        return ret;
    }

    public double[] solve(double[] input) {
        double[] transformFromIdentity = new double[input.length];
        System.arraycopy(input, 0, transformFromIdentity, 0, input.length);
        int n = coefficients.length;
        for (int i=0; i < n-1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                    int indexj = rowIndex[j];
                    transformFromIdentity[indexj]
                          -= coefficients[indexj][i] * transformFromIdentity[rowIndex[i]];
            }
        }

        double inverse[] = new double[n];
        inverse[n -1] =
              transformFromIdentity[rowIndex[n -1]]/coefficients[rowIndex[n -1]][n -1];
        for (int row = n -2; row>=0; --row)
        {
            inverse[row] = transformFromIdentity[rowIndex[row]];
            for (int k = row+1; k< n; ++k)
            {
                inverse[row] -= coefficients[rowIndex[row]][k]*inverse[k];
            }
            inverse[row] /= coefficients[rowIndex[row]][row];
        }
        return inverse;
    }

    public double[][] solveDoubleArray(double[][] transformFromIdentity) {
        int n = coefficients.length;
        for (int i=0; i < n-1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int COLumN = 0; COLumN < transformFromIdentity[0].length; ++COLumN) {
                    int indexj = rowIndex[j];
                    transformFromIdentity[indexj][COLumN]
                          -= coefficients[indexj][i] * transformFromIdentity[rowIndex[i]][COLumN];
                }
            }
        }

        double inverse[][] = new double[n][transformFromIdentity[0].length];
        for (int COLumN = 0; COLumN < transformFromIdentity[0].length; ++COLumN)
        {
            inverse[n -1][COLumN] =
                  transformFromIdentity[rowIndex[n -1]][COLumN]/coefficients[rowIndex[n -1]][n -1];
            for (int row = n -2; row>=0; --row)
            {
                inverse[row][COLumN] = transformFromIdentity[rowIndex[row]][COLumN];
                for (int k = row+1; k< n; ++k)
                {
                    inverse[row][COLumN] -= coefficients[rowIndex[row]][k]*inverse[k][COLumN];
                }
                inverse[row][COLumN] /= coefficients[rowIndex[row]][row];
            }
        }
        return inverse;
    }

    private double[][] populateTransformFromIdentity(double[][] values, int n) {
        double transformFromIdentity[][] = new double[n][n + values[0].length];
        for (int i = 0; i < n; ++i) {
            transformFromIdentity[i][i] = 1;
        }
        for (int Column = n; Column < transformFromIdentity[0].length; ++Column) {
            for (int row = 0; row < n; ++row) {
                transformFromIdentity[row][Column] = values[row][Column- n];
            }
        }
        return transformFromIdentity;
    }


    private double[][] populateValueHolder(double[][] values, int n) {
        double transformFromIdentity[][] = new double[n][values[0].length];
        for (int Column = 0; Column < transformFromIdentity[0].length; ++Column) {
            for (int row = 0; row < n; ++row) {
                transformFromIdentity[row][Column] = values[row][Column];
            }
        }
        return transformFromIdentity;
    }

    public  void gaussian(double coefficients[][])
    {
        int n = coefficients.length;
        rowIndex = new int[n];
        double colMax[] = new double[n];
 
        // Initialize the index
        for (int i = 0; i < n; ++i) {
            rowIndex[i] = i;
        }
 
        // Find the rescaling factors, one from each row
        for (int row = 0; row < n; ++row)
        {
            double max = 0;
            for (int col = 0; col<n; ++col)
            {
                double c0 = Math.abs(coefficients[row][col]);
                if (c0 > max) max = c0;
            }
            colMax[row] = max;
        }
 
        // Search the pivoting element from each column
        int k = 0;
        for (int column = 0; column < n-1; ++column)
        {
            double pivot_maxrow = 0;
            for (int i = column; i < n; ++i)
            {
                double pi0 = Math.abs(coefficients[rowIndex[i]][column]);
                pi0 /= colMax[rowIndex[i]];
                if (pi0 > pivot_maxrow)
                {
                    pivot_maxrow = pi0;
                    k = i;
                }
            }
 
            // Interchange rows according to the pivoting order
            int itmp = rowIndex[column];
            rowIndex[column] = rowIndex[k];
            rowIndex[k] = itmp;
            for (int i = column + 1; i < n; ++i)
            {
                double pj = coefficients[rowIndex[i]][column]/coefficients[rowIndex[column]][column];
 
                // Record pivoting ratios below the diagonal
                coefficients[rowIndex[i]][column] = pj;
 
                // Modify other elements accordingly (top triangular)
                for (int cols_right = column + 1; cols_right < n; ++cols_right) {
                    coefficients[rowIndex[i]][cols_right] -= pj * coefficients[rowIndex[column]][cols_right];
                }
            }
        }
    }
}