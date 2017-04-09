/**
 * 
 */
package math;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * NOT complete!!!
 * @author oshanker
 *
 */
public class ConjugateGradient2 {
	public static interface Function{
		double[] evaluate(double[] args);
	}
	static final double realmin = 2.2251e-308;
	
	/*
	 * Minimize a differentiable multivariate function. 
	 *
	 * Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, ... )
	 *
	 * where the starting point is given by "X" (D by 1), and the function named in
	 * the string "f", must return a function value and a vector of partial
	 * derivatives of f wrt X, the "length" gives the length of the run: if it is
	 * positive, it gives the maximum number of line searches, if negative its
	 * absolute gives the maximum allowed number of function evaluations. You can
	 * (optionally) give "length" a second component, which will indicate the
	 * reduction in function value to be expected in the first line-search (defaults
	 * to 1.0). The parameters P1, P2, P3, ... are passed on to the function f.
	 *
	 * The function returns when either its length is up, or if no further progress
	 * can be made (ie, we are at a (local) minimum, or so close that due to
	 * numerical problems, we cannot get any closer). NOTE: If the function
	 * terminates within a few iterations, it could be an indication that the
	 * function values and derivatives are not consistent (ie, there may be a bug in
	 * the implementation of your "f" function). The function returns the found
	 * solution "X", a vector of function values "fX" indicating the progress made
	 * and "i" the number of iterations (line searches or function evaluations,
	 * depending on the sign of "length") used.
	 *
	 * The Polack-Ribiere flavour of conjugate gradients is used to compute search
	 * directions, and a line search using quadratic and cubic polynomial
	 * approximations and the Wolfe-Powell stopping criteria is used together with
	 * the slope ratio method for guessing initial step sizes. Additionally a bunch
	 * of checks are made to make sure that exploration is taking place and that
	 * extrapolation will not be unboundedly large.
	 *
	 * See also: checkgrad 
	 *
	 * Copyright (C) 2001 - 2006 by Carl Edward Rasmussen (2006-09-08).
	 */
	
	public static int fmincg(Function f, double[] X, HashMap<String, String> options, ArrayList<Double> fX) {
		int length = 100;
		double red = 1; 

			// Read options
		if (options != null  ){
			if(options.containsKey("MaxIter")) {
			    length = Integer.parseInt(options.get("MaxIter"));
			}
			if(options.containsKey("red")) {
			    red = Double.parseDouble(options.get("red"));
			}
	    }
		/*
		 *  don't reevaluate within 0.1 of the limit of the current bracket
		 *  extrapolate maximum 3 times the current step-size
		 *  max 20 function evaluations per line search
		 *  maximum allowed slope ratio
		 *  SIG and RHO are the constants controlling the Wolfe-
		 *  Powell conditions. SIG is the maximum allowed absolute ratio between
		 *  previous and new slopes (derivatives in the search direction), thus setting
		 *  SIG to low (positive) values forces higher precision in the line-searches.
		 *  RHO is the minimum allowed fraction of the expected (from the slope at the
		 *  initial point in the linesearch). Constants must satisfy 0 < RHO < SIG < 1.
		 *  Tuning of SIG (depending on the nature of the function to be optimized) may
		 *  speed up the minimization; it is probably not worth playing much with RHO.

		 *  The code falls naturally into 3 parts, after the initial line search is
		 *  started in the direction of steepest descent. 1) we first enter a while loop
		 *  which uses point 1 (p1) and (p2) to compute an extrapolation (p3), until we
		 *  have extrapolated far enough (Wolfe-Powell conditions). 2) if necessary, we
		 *  enter the second loop which takes p2, p3 and p4 chooses the subinterval
		 *  containing a (local) minimum, and interpolates it, unil an acceptable point
		 *  is found (Wolfe-Powell conditions). Note, that points are always maintained
		 *  in order p0 <= p1 <= p2 < p3 < p4. 3) compute a new search direction using
		 *  conjugate gradients (Polack-Ribiere flavour), or revert to steepest if there
		 *  was a problem in the previous line-search. Return the best value so far, if
		 *  two consecutive line-searches fail, or whenever we run out of function
		 *  evaluations or line-searches. During extrapolation, the "f" function may fail
		 *  either with an error or returning Nan or Inf, and minimize should handle this
		 *  gracefully.
		 */
		double SIG = 0.1;       // RHO and SIG are the constants in the Wolfe-Powell conditions
		double RHO = SIG/2;     // a bunch of constants for line searches
		double INT = 0.1;    // don't reevaluate within 0.1 of the limit of the current bracket
		double EXT = 3.0;                    // extrapolate maximum 3 times the current bracket
		int MAX = 20;                         // max 20 function evaluations per line search
		double RATIO = 10;                                      // maximum allowed slope ratio
		String S="Iteration";

		int i = 0;                                            // zero the run length counter
		boolean ls_failed = false;                             // no previous line search has failed
		//fX = [];
		double[] f0df0_new = f.evaluate(X); 
		// get function value and gradient
		double f0_new = f0df0_new[0];
		double[] df0_new = new double[f0df0_new.length-1];
		System.arraycopy(f0df0_new, 1, df0_new, 0, df0_new.length);
		if(fX != null) { fX.add(f0_new);}
		if(length<0){i++;}                                            // count epochs?!
		double[] s = new double[df0_new.length];
		System.arraycopy(df0_new, 0, s, 0, df0_new.length);
		// search direction is steepest
		for (int j = 0; j < s.length; j++) {
			s[i] = -s[i];
		}
		double d0_new = -arrayProduct(s, s);
		// this is the slope
		double x3_new = red/(1-d0_new);                                  // initial step is red/(|s|+1)
		while( i < Math.abs(length)){      // while not finished
			if(length>0){ i++ ;   }        // count iterations?!
			double[] X0 = new double[X.length];
			System.arraycopy(X, 0, X0, 0, X.length);

			double F0_new = f0_new; 
			double[] dF0_new = new double[df0_new.length];
			System.arraycopy(df0_new, 0, dF0_new, 0, df0_new.length);
			// make a copy of current values (line 79)
			int M;
			if (length>0){ M = MAX;} else{ M = Math.min(MAX, -length-i);}
			boolean success = false; 
			while (true){ //keep extrapolating as long as necessary 82
			    //x2 = 0; f2 = f0; d2 = d0; f3 = f0; df3 = df0; TODO 83
				success = false; 
				double A;
				double B;
				while (!success && (M > 0)) {
					M = M - 1; if(length<0){i++;}   // count epochs?!
				}//end while !success && (M > 0)
			}//end  while true // end extrapolation
			
			/*

  while (abs(d3) > -SIG*d0 || f3 > f0+x3*RHO*d0) && M > 0  % keep interpolating
    if d3 > 0 || f3 > f0+x3*RHO*d0                         % choose subinterval
      x4 = x3; f4 = f3; d4 = d3;                      % move point 3 to point 4
    else
      x2 = x3; f2 = f3; d2 = d3;                      % move point 3 to point 2
    end
    if f4 > f0           
      x3 = x2-(0.5*d2*(x4-x2)^2)/(f4-f2-d2*(x4-x2));  % quadratic interpolation
    else
      A = 6*(f2-f4)/(x4-x2)+3*(d4+d2);                    % cubic interpolation
      B = 3*(f4-f2)-(2*d2+d4)*(x4-x2);
      x3 = x2+(sqrt(B*B-A*d2*(x4-x2)^2)-B)/A;        % num. error possible, ok!
    end
    if isnan(x3) || isinf(x3)
      x3 = (x2+x4)/2;               % if we had a numerical problem then bisect
    end
    x3 = max(min(x3, x4-INT*(x4-x2)),x2+INT*(x4-x2));  % don't accept too close
    [f3 df3] = feval(f, X+x3*s, varargin{:});
    if f3 < F0, X0 = X+x3*s; F0 = f3; dF0 = df3; end         % keep best values
    M = M - 1; i = i + (length<0);                             % count epochs?!
    d3 = df3'*s;                                                    % new slope
  end                                                       % end interpolation

  if abs(d3) < -SIG*d0 && f3 < f0+x3*RHO*d0          % if line search succeeded
    X = X+x3*s; f0 = f3; fX = [fX' f0]';                     % update variables
    fprintf('%s %6i;  Value %4.6e\r', S, i, f0);
    s = (df3'*df3-df0'*df3)/(df0'*df0)*s - df3;   % Polack-Ribiere CG direction
    df0 = df3;                                               % swap derivatives
    d3 = d0; d0 = df0'*s;
    if d0 > 0                                      % new slope must be negative
      s = -df0; d0 = -s'*s;                  % otherwise use steepest direction
    end
    x3 = x3 * min(RATIO, d3/(d0-realmin));          % slope ratio but max RATIO
    ls_failed = 0;                              % this line search did not fail
  else
    X = X0; f0 = F0; df0 = dF0;                     % restore best point so far
    if ls_failed || i > abs(length)         % line search failed twice in a row
      break;                             % or we ran out of time, so we give up
    end
    s = -df0; d0 = -s'*s;                                        % try steepest
    x3 = 1/(1-d0);                     
    ls_failed = 1;                                    % this line search failed
  end

			
		*/
		}//end while i  < Math.abs(length)
		//fprintf('\n');
		return i;
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
			public double[] evaluate(double[] args) {
				double[] ret = new double[args.length+1];
				for (int i = 0; i < args.length; i++) {
					ret[0] += (args[i]-1)*(args[i]-1);
					ret[i+1] = 2*(args[i]-1);
				}
				return ret;
			}
		};
		double[] X = {1,2};
		ArrayList<Double> fX = new ArrayList<>();
		HashMap<String, String> options = null;
		int ans = fmincg(f, X , options , fX);
		System.out.println("i " + ans + " " + Arrays.toString(X));
		System.out.println(fX);
	}

}
