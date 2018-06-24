df <- read.csv("../oldriemann/out/gzetaE12/gzeta6.csv", header = FALSE);
i = 10;
	x<-as.matrix(df[i]);
	x<-as.vector(x);
b=seq(-10,10,0.5)
yy = range(x)
b=c(yy[1],b,yy[2])
h=hist(x,breaks=b,xlim=c(-10,10),xlab='Z',
main = paste("Histogram of Z for 3" , expression(pi), "/2"))
