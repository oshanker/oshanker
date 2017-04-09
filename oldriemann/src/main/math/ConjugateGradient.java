/**
 * 
 */
package math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * @author oshanker
 *
 */
public class ConjugateGradient {
	public static interface Function{
		double[] eval(double[] args);
	}
	static final double realmin = 2.2251e-308;
	
	//function [X, fX, i] = fmincg(f, X, options)
			// Minimize a continuous differentialble multivariate function. Starting point
			// is given by "X" (D by 1), and the function named in the string "f", must
			// return a function value and a vector of partial derivatives. The Polack-
			// Ribiere flavour of conjugate gradients is used to compute search directions,
			// and a line search using quadratic and cubic polynomial approximations and the
			// Wolfe-Powell stopping criteria is used together with the slope ratio method
			// for guessing initial step sizes. Additionally a bunch of checks are made to
			// make sure that exploration is taking place and that extrapolation will not
			// be unboundedly large. The "length" gives the length of the run: if it is
			// positive, it gives the maximum number of line searches, if negative its
			// absolute gives the maximum allowed number of function evaluations. You can
			// (optionally) give "length" a second component, which will indicate the
			// reduction in function value to be expected in the first line-search (defaults
			// to 1.0). The function returns when either its length is up, or if no further
			// progress can be made (ie, we are at a minimum, or so close that due to
			// numerical problems, we cannot get any closer). If the function terminates
			// within a few iterations, it could be an indication that the function value
			// and derivatives are not consistent (ie, there may be a bug in the
			// implementation of your "f" function). The function returns the found
			// solution "X", a vector of function values "fX" indicating the progress made
			// and "i" the number of iterations (line searches or function evaluations,
			// depending on the sign of "length") used.
			//
			// Usage: [X, fX, i] = fmincg(f, X, options, P1, P2, P3, P4, P5)
			//
			// See also: checkgrad 
			//
			// Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13
			//
			//
			// (C) Copyright 1999, 2000 & 2001, Carl Edward Rasmussen
			// 
			// Permission is granted for anyone to copy, use, or modify these
			// programs and accompanying documents for purposes of research or
			// education, provided this copyright notice is retained, and note is
			// made of any changes that have been made.
			// 
			// These programs and documents are distributed without any warranty,
			// express or implied.  As the programs were written for research
			// purposes only, they have not been tested to the degree that would be
			// advisable in any important application.  All use of these programs is
			// entirely at the user's own risk.
			//
			// [ml-class] Changes Made:
			// 1) Function name and argument specifications
			// 2) Output display
			//
	
	public static double[] fmincg(Function f, double[] X, HashMap<String, String> options) {
		int length = 100;

			// Read options
		if (options != null && options.containsKey("MaxIter")){
			 length = Integer.parseInt(options.get("MaxIter"));
	    }
		double RHO = 0.01;                            // a bunch of constants for line searches
		double SIG = 0.5;       // RHO and SIG are the constants in the Wolfe-Powell conditions
		double INT = 0.1;    // don't reevaluate within 0.1 of the limit of the current bracket
		double EXT = 3.0;                    // extrapolate maximum 3 times the current bracket
		int MAX = 20;                         // max 20 function evaluations per line search
		double RATIO = 100;                                      // maximum allowed slope ratio
		double red = 1; 
		String S="Iteration";
		ArrayList<Double> fX = new ArrayList<>();

		int i = 0;                                            // zero the run length counter
		boolean ls_failed = false;                             // no previous line search has failed
		//fX = [];
		double[] f1df1 = f.eval(X); 
		// get function value and gradient
		double f1 = f1df1[0];
		double[] df1 = new double[f1df1.length-1];
		System.arraycopy(f1df1, 1, df1, 0, df1.length);
		if(length<0){i++;}                                            // count epochs?!
		double[] s = new double[df1.length];
		System.arraycopy(df1, 0, s, 0, df1.length);
		// search direction is steepest
		for (int j = 0; j < s.length; j++) {
			s[i] = -s[i];
		}
		double d1 = -arrayProduct(s, s);
		// this is the slope
		double z1 = red/(1-d1);                                  // initial step is red/(|s|+1)
		while( i < Math.abs(length)){      // while not finished
			if(length>0){ i++ ;   }        // count iterations?!
			double[] X0 = new double[X.length];
			System.arraycopy(X, 0, X0, 0, X.length);

			double f0 = f1; 
			double[] df0 = new double[df1.length];
			System.arraycopy(df1, 0, df0, 0, df1.length);
			// make a copy of current values

			for (int j = 0; j < s.length; j++) {
				X[j] += z1 *s[j];
			}
			// begin line search
			double[] f2df2 = f.eval(X); 
			// get function value and gradient
			double f2 = f2df2[0];
			double[] df2 = new double[f2df2.length-1];
			System.arraycopy(f2df2, 1, df2, 0, df2.length);
			if(length<0){i++;}  // count epochs?!
			double d2 = arrayProduct(s, df2);
			double f3 = f1; double d3 = d1; double z3 = -z1;  // initialize point 3 equal to point 1
			int M;
			if (length>0){ M = MAX;} else{ M = Math.min(MAX, -length-i);}
			boolean success = false; 
			double limit = -1;                     // initialize quanteties
			while (true){
				double z2;
				double A;
				double B;
				while (((f2 > f1+z1*RHO*d1) | (d2 > -SIG*d1)) & (M > 0) ){
					limit = z1;                                         // tighten the bracket
					if (f2 > f1) {
						z2 = z3 - (0.5*d3*z3*z3)/(d3*z3+f2-f3);                 // quadratic fit
					} else {
						A = 6*(f2-f3)/z3+3*(d2+d3);                                 // cubic fit
						B = 3*(f3-f2)-z3*(d3+2*d2);
						z2 = (Math.sqrt(B*B-A*d2*z3*z3)-B)/A;       // numerical error possible - ok!
					} //end
					//Double.
					if( Double.isFinite(z2) ) {
						z2 = z3/2;                  // if we had a numerical problem then bisect
					}
					z2 = Math.max(Math.min(z2, INT*z3),(1-INT)*z3);  // don't accept too close to limits
					z1 = z1 + z2;                                           // update the step
					for (int j = 0; j < X.length; j++) {
						X[j] = X[j] + z2*s[j];
					}
					f2df2 = f.eval(X); 
					// get function value and gradient
					f2 = f2df2[0];
					//df2 = new double[f2df2.length-1];
					System.arraycopy(f2df2, 1, df2, 0, df2.length);

					//[f2 df2] = eval(argstr);
					// fprintf('quadratic //s //4i | Cost: //4.6e\n', S, i, f2);
					M = M - 1; if(length<0){i++;}   // count epochs?!
					d2 = arrayProduct(s, df2);
					z3 = z3-z2;                    // z3 is now relative to the location of z2
				}//end while (((f2 > 
				if (f2 > f1+z1*RHO*d1 | d2 > -SIG*d1) {
					break;                                                // this is a failure
				} else if (d2 > SIG*d1) {
					success = true; break;                                             // success
				} else if (M == 0) {
					break;                                                          // failure
				} //end
				A = 6*(f2-f3)/z3+3*(d2+d3);                      // make cubic extrapolation
				B = 3*(f3-f2)-z3*(d3+2*d2);
				z2 = -d2*z3*z3/(B+Math.sqrt(B*B-A*d2*z3*z3));        // num. error possible - ok!
				if (Double.isFinite(z2) | z2 < 0 ) {  // num prob or wrong sign?
					if (limit < -0.5 ) {                             // if we have no upper limit
						z2 = z1 * (EXT-1);                 // the extrapolate the maximum amount
					} else {
						z2 = (limit-z1)/2;                                   // otherwise bisect
					} //end
				} else if( (limit > -0.5) & (z2+z1 > limit) ){  // extraplation beyond max?
					z2 = (limit-z1)/2;                                               // bisect
				} else if( (limit < -0.5) & (z2+z1 > z1*EXT) ){      // extrapolation beyond limit
					z2 = z1*(EXT-1.0);                           // set to extrapolation limit
				} else if( z2 < -z3*INT ){ 
					z2 = -z3*INT;
				} else if( (limit > -0.5) & (z2 < (limit-z1)*(1.0-INT)) ){  // too close to limit?
					z2 = (limit-z1)*(1.0-INT);
				}//end
				f3 = f2; d3 = d2; z3 = -z2;                  // set point 3 equal to point 2
				z1 = z1 + z2; //X = X + z2*s;   // update current estimates
				for (int j = 0; j < X.length; j++) {
					X[j] = X[j] + z2*s[j];
				}
				f2df2 = f.eval(X); // get function value and gradient
				f2 = f2df2[0];
				//df2 = new double[f2df2.length-1];
				System.arraycopy(f2df2, 1, df2, 0, df2.length);
				//fprintf('cubic //s //4i | Cost: //4.6e\n', S, i, f2);
				M = M - 1; if(length<0){i++;}   // count epochs?!
				d2 = arrayProduct(s, df2);
			}//end  while true // end of line search

			if (success){                 // if line search succeeded
				f1 = f2; fX.add(f1);
				double tmp = (arrayProduct(df2, df2)-arrayProduct(df1, df2))/arrayProduct(df1, df1);
				for (int j = 0; j < df2.length; j++) {
					s[j] = tmp*s[j] - df2[j]; // Polack-Ribiere direction
				}
				/*
			    fprintf('//s //4i | Cost: //4.6e\r', S, i, f1);
			    s = (df2'*df2-df1'*df2)/(df1'*df1)*s - df2;      // Polack-Ribiere direction
				 */
				swap(df1,df2);    // swap derivatives
				d2 = arrayProduct(df1,s);
				if (d2 > 0 ) {                // new slope must be negative
					//s = -df1;          // otherwise use steepest direction
					for (int j = 0; j < df1.length; j++) {
						s[j] =  - df1[j]; // otherwise use steepest direction
					}
					d2 = -arrayProduct(s,s);    
				}// end
				z1 = z1 * Math.min(RATIO, d1/(d2-realmin));          // slope ratio but max RATIO
				d1 = d2;
				ls_failed = false;            // this line search did not fail
			} else {
				//X = X0; f1 = f0; df1 = df0;  // restore point from before failed line search
				f1 = f0;
				System.arraycopy(X0, 0, X, 0, X.length);
				System.arraycopy(df0, 0, df1, 0,df1.length);
				
			    if (ls_failed | i > Math.abs(length) ){         // line search failed twice in a row
			      break;                             // or we ran out of time, so we give up
			    } //end
			    swap(df1,df2);    // swap derivatives
				for (int j = 0; j < df1.length; j++) { // try steepest
					s[j] =  - df1[j]; // otherwise use steepest direction
				}
				d1 = -arrayProduct(s,s);    
			    z1 = 1/(1-d1);                     
			    ls_failed = true;    // this line search failed
			}//end
		}//end while i  < Math.abs(length)
		//fprintf('\n');
		return X;
    }

	private static double arrayProduct(double[] s, double[] df2) {
		double d2 = 0;
		// d2 = df2'*s;
		for (int j = 0; j < s.length; j++) {
			d2 += df2[j] *s[j];
		}
		return d2;
	}

	private static void swap(double[] s, double[] df2) {
		for (int j = 0; j < s.length; j++) {
		    double tmp = s[j]; s[j] = df2[j]; df2[j] = tmp;    // swap derivatives
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Function f = new Function() {
			public double[] eval(double[] args) {
				double[] ret = new double[args.length+1];
				for (int i = 0; i < args.length; i++) {
					ret[0] += (args[i]-1)*(args[i]-1);
					ret[i+1] = 2*(args[i]-1);
				}
				return ret;
			}
		};
		double[] X = {1,2};
		double[] ans = fmincg(f, X , null);
		System.out.println(Arrays.toString(ans));
	}

}
