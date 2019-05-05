basedir<-"../oldriemann/out/gzetaE28"
df <- read.csv(paste0(basedir,"/gzeta12.csv"), header = FALSE);
outtable<-matrix(nrow =  12, ncol = length(df), byrow = TRUE)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(quantile(x, probs = seq(0, 1, 0.1)))
	summary <- c((i-1), summary)
	outtable[1:length(summary), i]<-summary
}
write(outtable,file=paste0(basedir,"/percentile.txt"),nc=dim(outtable)[1])
