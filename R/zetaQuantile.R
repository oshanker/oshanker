basedir<-"../oldriemann/out/gzetaE28"
df <- read.csv(paste0(basedir,"/gzeta12.csv"), header = FALSE);
outtable<-matrix(nrow =  length(df), ncol = 12, byrow = TRUE)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(quantile(x, probs = seq(0, 1, 0.1)))
	summary <- c((i-1), summary)
	outtable[i, 1:length(summary)]<-summary
}
write.csv(outtable,file=paste0(basedir,"/percentile.csv"))
mytest=read.csv(paste0(basedir,"/percentile.csv"))
print(mytest)
