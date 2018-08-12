df <- read.csv("../oldriemann/out/gzetaE12/gzeta6.csv", header = FALSE);
library(moments)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	s = skewness(x)
	k = kurtosis(x)
	summary = c(i-1, s, k, sd(x))
    conv = format( summary, scientific = FALSE, drop0trailing = TRUE, 
              digits = 4, width =8)
    print(conv, quote = FALSE)
}
