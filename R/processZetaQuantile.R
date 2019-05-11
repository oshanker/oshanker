basedir<-"../oldriemann/out/gzetaE12"
#basedir<-"../oldriemann/out/gzetaE28"
intable <- read.csv(paste0(basedir,"/percentile_calc.csv"),header=TRUE);
#intable <- read.csv(paste0(basedir,"/percentileE28.csv"),header=TRUE);
outfile<-"/quantileFitted_calc6.csv"
#outfile<-"/quantileFittedE28.csv"

colcount = dim(intable)[2]
print(paste("colcount",colcount))
fitted<-matrix(nrow=3,ncol=colcount)
meanValues<- c(1:colcount)
slopes<- c(1:colcount)
	x<-as.matrix(intable[colcount]);
	x<-as.vector(x);

for (i in 2:(colcount-1)) {
	y<-as.matrix(intable[i]);
	y<-as.vector(y);
#	print(x)
#	print(paste('is',i))
	meanValues[i]=mean(y) 
	slopes[i]=(y[3]-y[5])/2
	
	lm.fit <- lm(y~x)
	fitted[1:2,i]<-coef(lm.fit)
	z<-summary(lm.fit)
	fitted[3,i]<-z$r.squared
}
    
conv = format( meanValues, scientific = FALSE, drop0trailing = TRUE, digits = 4, width =8)
print(conv, quote = FALSE)
print(slopes)

write.csv(fitted,file=paste0(basedir, outfile), row.names = FALSE)
