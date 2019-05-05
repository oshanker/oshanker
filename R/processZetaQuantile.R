basedir<-"../oldriemann/out/gzetaE28"
outtable <- read.table(paste0(basedir,"/percentile.txt"), header = FALSE);
meanValues<- c(1:12)
slopes<- c(1:12)
for (i in 2:12) {
	x<-as.matrix(outtable[i]);
	x<-as.vector(x);
	meanValues[i]=mean(x) 
	slopes[i]=(x[1]-x[13])/4
}
print(meanValues)
print(slopes)