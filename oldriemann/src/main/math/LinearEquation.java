package math;

import java.util.Arrays;

public class LinearEquation
{
    double [][] coefficients;
    double [][] values;
    int index[];

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

    void printMatrix(double[][] inverted_mat) {
        for (int i = 0; i < inverted_mat.length; ++i)
        {
            for (int j = 0; j < inverted_mat[i].length; ++j)
            {
                System.out.print(inverted_mat[i][j] + "  ");
            }
            System.out.println();
        }
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

    public LinearEquation(double[][] coefficients) {
        this.coefficients = coefficients;
        gaussian(coefficients);
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

    public double[] solve(double[] transformFromIdentity) {
        int n = coefficients.length;
        for (int i=0; i < n-1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                    int indexj = index[j];
                    transformFromIdentity[indexj]
                          -= coefficients[indexj][i] * transformFromIdentity[index[i]];
            }
        }

        double inverse[] = new double[n];
        inverse[n -1] =
              transformFromIdentity[index[n -1]]/coefficients[index[n -1]][n -1];
        for (int row = n -2; row>=0; --row)
        {
            inverse[row] = transformFromIdentity[index[row]];
            for (int k = row+1; k< n; ++k)
            {
                inverse[row] -= coefficients[index[row]][k]*inverse[k];
            }
            inverse[row] /= coefficients[index[row]][row];
        }
        return inverse;
    }

    public double[][] solveDoubleArray(double[][] transformFromIdentity) {
        int n = coefficients.length;
        for (int i=0; i < n-1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int COLumN = 0; COLumN < transformFromIdentity[0].length; ++COLumN) {
                    int indexj = index[j];
                    transformFromIdentity[indexj][COLumN]
                          -= coefficients[indexj][i] * transformFromIdentity[index[i]][COLumN];
                }
            }
        }

        double inverse[][] = new double[n][transformFromIdentity[0].length];
        for (int COLumN = 0; COLumN < transformFromIdentity[0].length; ++COLumN)
        {
            inverse[n -1][COLumN] =
                  transformFromIdentity[index[n -1]][COLumN]/coefficients[index[n -1]][n -1];
            for (int row = n -2; row>=0; --row)
            {
                inverse[row][COLumN] = transformFromIdentity[index[row]][COLumN];
                for (int k = row+1; k< n; ++k)
                {
                    inverse[row][COLumN] -= coefficients[index[row]][k]*inverse[k][COLumN];
                }
                inverse[row][COLumN] /= coefficients[index[row]][row];
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
        index = new int[n];
        double c[] = new double[n];
 
        // Initialize the index
        for (int i=0; i<n; ++i) {
            index[i] = i;
        }
 
        // Find the rescaling factors, one from each row
        for (int i=0; i<n; ++i) 
        {
            double max = 0;
            for (int j=0; j<n; ++j) 
            {
                double c0 = Math.abs(coefficients[i][j]);
                if (c0 > max) max = c0;
            }
            c[i] = max;
        }
 
        // Search the pivoting element from each column
        int k = 0;
        for (int j=0; j<n-1; ++j) 
        {
            double pivot_maxrow = 0;
            for (int i = j; i < n; ++i)
            {
                double pi0 = Math.abs(coefficients[index[i]][j]);
                pi0 /= c[index[i]];
                if (pi0 > pivot_maxrow)
                {
                    pivot_maxrow = pi0;
                    k = i;
                }
            }
 
            // Interchange rows according to the pivoting order
            int itmp = index[j];
            index[j] = index[k];
            index[k] = itmp;
            for (int i=j+1; i<n; ++i) 	
            {
                double pj = coefficients[index[i]][j]/coefficients[index[j]][j];
 
                // Record pivoting ratios below the diagonal
                coefficients[index[i]][j] = pj;
 
                // Modify other elements accordingly
                for (int l=j+1; l<n; ++l) {
                    coefficients[index[i]][l] -= pj * coefficients[index[j]][l];
                }
            }
        }
    }
}