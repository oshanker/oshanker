df <- read.csv("../oldriemann/out/gzetaE28/gzeta6.csv", header = FALSE);
library(moments)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	s = skewness(x)
	k = kurtosis(x)
	print(c(i-1, s, k))
}
