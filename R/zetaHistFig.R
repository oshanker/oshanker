df <- read.csv("../oldriemann/out/gzetaE12/gzeta6.csv", header = FALSE);
i = 10;
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	x<-log(abs(x));
library(moments)
#qqnorm(x, main = paste("Normal Q-Q Plot for ln(|Z(t)|) at ", expression(phi),
#   " = 3", expression(pi), "/2"),
#       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); qqline(x, col = 2)
	s = skewness(x)
	k = kurtosis(x)
	print(c( s, k))
       
b=seq(-10,10,0.5)
yy = range(x)
#b=c(yy[1],b,yy[2])
#h=hist(x,breaks=b,xlim=c(-10,10),xlab='Z',
# main = paste("Histogram of Z for 3" , expression(pi), "/2"))
