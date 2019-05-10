basedir<-"../oldriemann/out/gzetaE28"
intable <- read.csv(paste0(basedir,"/percentileE28.csv"),header=TRUE);
meanValues<- c(1:13)
slopes<- c(1:13)
for (i in 2:13) {
	x<-as.matrix(intable[i]);
	x<-as.vector(x);
#	print(x)
#	print(paste('is',i))
	meanValues[i]=mean(x) 
	slopes[i]=(x[3]-x[5])/2
}
    conv = format( meanValues, scientific = FALSE, drop0trailing = TRUE, digits = 4, width =8)
    print(conv, quote = FALSE)
#print(meanValues)
print(slopes)

 outtable<-matrix(nrow=2,ncol=length(meanValues))
 outtable[1,1:dim(outtable)[2]]<-meanValues
 outtable[2,1:dim(outtable)[2]]<-slopes
 write.csv(outtable,file=paste0(basedir,"/quantileSlopesE28.csv"), row.names = FALSE)
