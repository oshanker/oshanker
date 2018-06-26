df <- read.csv("../oldriemann/out/gzetaE12/gzeta.csv", header = FALSE);
library(moments)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	s = skewness(x)
	k = kurtosis(x)
	print(c(i, s, k))
}
