df <- read.csv("../oldriemann/out/gzetaE12/gzeta6.csv", header = FALSE);
i = 10;
	x<-as.matrix(df[i]);
	x<-as.vector(x);
qqnorm(x, main = paste("Normal Q-Q Plot for Z(t) at ", expression(phi),
   " = 3", expression(pi), "/2"),
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); qqline(x, col = 2)
       
b=seq(-10,10,0.5)
yy = range(x)
#b=c(yy[1],b,yy[2])
#h=hist(x,breaks=b,xlim=c(-10,10),xlab='Z',
# main = paste("Histogram of Z for 3" , expression(pi), "/2"))
