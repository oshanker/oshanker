basedir<-"../oldriemann/out/gzetaE28"
intable <- read.csv(paste0(basedir,"/percentile.csv"));
outtable<-matrix(nrow =  2, ncol = 12, byrow = TRUE)
meanValues<- c(1:12)
slopes<- c(1:12)
for (i in 2:12) {
	x<-as.matrix(intable[i]);
	x<-as.vector(x);
	meanValues[i]=mean(x) 
	slopes[i]=(x[4]-x[10])/2.828
}
outtable[1,1:length(meanValues)]<-meanValues
outtable[2,1:length(slopes)]<-slopes
write.csv(outtable,file=paste0(basedir,"/quantileSlopes.csv"))
mytest=read.csv(paste0(basedir,"/quantileSlopes.csv"))
print(mytest)
