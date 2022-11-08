package math;

import java.util.Scanner;
 
public class LinearEquation 
{
    int n = 5;
    double [][] coefficients = new double[n][n];
    double [][] values = new double[n][1];

    public static void main(String args[]) {
        LinearEquation linearEquation = new LinearEquation();
        linearEquation.runInvert();
    }

    public void runInvert() {
        initEquations();
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
        double inverted_mat[][] = invert(coefficients, values);
        System.out.println("The inverse is: ");
        for (int i=0; i< inverted_mat.length; ++i)
        {
            for (int j=0; j < inverted_mat[i].length; ++j)
            {
                System.out.print(inverted_mat[i][j] + "  ");
            }
            System.out.println();
        }
        //Multiplication of mat inverse and constants
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

    private void initEquations() {
        for(int row = 0; row < n; row++)
        {
            for(int col = 0; col <= row; col ++)
            {
                coefficients[row][col ] = 1;
            }
            values[row][0] = row + 1.5;
        }
        coefficients[0][n-2 ] += 80;
        values[0][0] += 80;
        coefficients[1][n-2 ] += 70;
        values[1][0] += 70;
        coefficients[n-1][n-1 ] += 20;
        values[n-1][0] += 41;
    }

    public static double[][] invert(double coefficients[][], double [][] values)
    {
        int n = coefficients.length;
        double transformFromIdentity[][] = new double[n][n+values[0].length];
        int index[] = new int[n];
        for (int i=0; i < n; ++i) {
            transformFromIdentity[i][i] = 1;
        }
        for (int Column=n; Column < transformFromIdentity[0].length; ++Column) {
            for (int row=0; row < n; ++row) {
                transformFromIdentity[row][Column] = values[row][Column-n];
            }
        }

        // Transform the matrix into an upper triangle
        gaussian(coefficients, index);
 
        // Update the matrix b[i][j] with the ratios stored
        for (int i=0; i<n-1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                for (int COLumN = 0; COLumN < transformFromIdentity[0].length; ++COLumN) {
                    int indexj = index[j];
                    transformFromIdentity[indexj][COLumN]
                          -= coefficients[indexj][i] * transformFromIdentity[index[i]][COLumN];
                }
            }
        }

        double inverse[][] = new double[n][n+values[0].length];
        // Perform backward substitutions
        for (int COLumN  = 0; COLumN < transformFromIdentity[0].length; ++COLumN)
        {
            inverse[n-1][COLumN] =
                  transformFromIdentity[index[n-1]][COLumN]/coefficients[index[n-1]][n-1];
            for (int row = n-2; row>=0; --row)
            {
                inverse[row][COLumN] = transformFromIdentity[index[row]][COLumN];
                for (int k=row+1; k<n; ++k)
                {
                    inverse[row][COLumN] -= coefficients[index[row]][k]*inverse[k][COLumN];
                }
                inverse[row][COLumN] /= coefficients[index[row]][row];
            }
        }
        return inverse;
    }
 

    public static void gaussian(double coefficients[][], int index[])
    {
        int n = index.length;
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