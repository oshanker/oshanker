basedir<-"../oldriemann/out/gzetaE12"
intable <- read.csv(paste0(basedir,"/percentile_calc.csv"),header=TRUE);
fitted<-matrix(nrow=3,ncol=13)
meanValues<- c(1:13)
slopes<- c(1:13)
	x<-as.matrix(intable[13]);
	x<-as.vector(x);

for (i in 2:12) {
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

 outtable<-matrix(nrow=2,ncol=length(meanValues))
 outtable[1,1:dim(outtable)[2]]<-meanValues
 outtable[2,1:dim(outtable)[2]]<-slopes
 write.csv(outtable,file=paste0(basedir,"/quantileSlopes_calc6.csv"), row.names = FALSE)
write.csv(fitted,file=paste0(basedir,"/quantileFitted_calc6.csv"), row.names = FALSE)
