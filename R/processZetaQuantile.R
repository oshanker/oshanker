#input form zetaQuantile.R
basedir<-"../oldriemann/out/gzetaE12"
#basedir<-"../oldriemann/out/gzetaE28"
intable <- read.csv(paste0(basedir,"/percentile_calc.csv"),header=TRUE);
#intable <- read.csv(paste0(basedir,"/percentileE28.csv"),header=TRUE);
outfile<-"/quantileFitted_calc6.csv"
#outfile<-"/quantileFittedE28.csv"

rowcount = dim(intable)[1]
colcount = dim(intable)[2]
print(paste("colcount",colcount))
fitted<-matrix(nrow=3,ncol=colcount-3)
meanValues<- c(1:colcount)
slopes<- c(1:colcount)
x<-as.matrix(intable[1]);
x<-as.vector(x);
x=x[2:rowcount]
x = cos(2*x*pi/rowcount)
print(x)
for (i in 3:colcount-1) {
	y<-as.matrix(intable[i]);
	y<-as.vector(y);
	y=y[2:rowcount]
	print( y)
#	print(x)
#	print(paste('is',i))

	lm.fit <- lm(y~x)
	fitted[1:2,i-3]<-coef(lm.fit)
	z<-summary(lm.fit)
	print(z)
	fitted[3,i-3]<-z$r.squared
}
    
print(fitted)
write.csv(fitted,file=paste0(basedir, outfile), row.names = FALSE)
