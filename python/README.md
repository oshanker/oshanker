# oshanker

We present here the value distribution for the Riemann zeta function and for RMT. We present the distribution for Generalized Gram point with $\phi$ = 90 or 270 degrees. More precisely, we are talking about the function $A$ discussed in [Random Matrix Theory explanation for Riemann Zeta Value Distribution Symmetry
](https://www.researchgate.net/publication/373331662_Random_Matrix_Theory_explanation_for_Riemann_Zeta_Value_Distribution_Symmetry). But they are close enough to each other. The distribution is the sum of a gaussian normal distribution and an exponential distribution (exp(- $lambda$ |x|)). Note that the exponential distribution has a discontinuous derivative at x = 0. 

The definitions, in detail, are


   gauss_norm = 1/(math.sqrt(2*math.pi))


def gauss(x, sigma ):

    A = gauss_norm/(sigma)
    
    return A*np.exp(-(x)**2/(2*sigma**2))

def exp_pdf(x, lam):

    xx = np.abs(np.copy(x))
    
    return  np.exp(-lam * xx) * lam/2
    
def exp_gauss(x, p, lam, sigma):

    ret = p * exp_pdf(x, lam) + (1-p) * gauss(x, sigma )
    
    return ret
    
See [characteristic.py](https://github.com/oshanker/oshanker/blob/master/python/cue/characteristic.py).

N  | T      | p    | $lambda$ | $sigma$ | source | 
---| ---    | ---  | ---      |---      |--- |
NA | 1.6E9  | 0.70 | 0.87     | 5.96    | Dirichlet | 
NA | 1.0E12 | 0.66 | 0.946    | 6.09    | Riemann | 
95 | 1.1E42 | 0.29 | 3.65     | 1.51    | CUE | 
150 | 8.75E65 | 0.285 | 4.23     | 1.46    | CUE | 


param [0.43206767, 2.04694732, 3.56966437]
vars [77.52212123 77.42877897 77.27782274 77.07953997 76.84744331 76.59734974
 76.34630274 76.11141075 75.90868126 75.75192995 75.65183918 75.61522995
 75.64459713 75.73793939 75.88889562 76.08717839 76.31927505 76.56936862
 76.82041562 77.05530761 77.2580371  77.41478841 77.51487918 77.55148841]

[Exponential distribution - Wikipedia](https://en.wikipedia.org/wiki/Exponential_distribution)

[algorithm - Exponential Distribution](https://stackoverflow.com/questions/2106503/pseudorandom-number-generator-exponential-distribution)

[Generate random numbers according to distributions](https://stackoverflow.com/questions/3510475/generate-random-numbers-according-to-distributions)

[fit-a-gaussian-curve](https://stackoverflow.com/questions/44480137/how-can-i-fit-a-gaussian-curve-in-python)

[fit-mixture-of-two-gaussian](https://stackoverflow.com/questions/35990467/fit-mixture-of-two-gaussian-normal-distributions-to-a-histogram-from-one-set-of)

[scipy.optimize.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html)

[area under normal curve](https://discovery.cs.illinois.edu/learn/Simulation-and-Distributions/Normal-Distribution/)

curl 'https://bootstrap.pypa.io/get-pip.py' -o get-pip.py

python get-pip.py --user

pip install pandas

python plot_quantiles.py 


https://machinelearningmastery.com/start-here/#deep_learning_time_series
