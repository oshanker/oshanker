basedir<-"../oldriemann/out/gzetaE28"
#lnfactor<-sqrt(log(1.0E12/(2*pi)))
df <- read.csv(paste0(basedir,"/gzeta6.csv"), header = FALSE);
outtable<-matrix(nrow =  length(df), ncol = 13, byrow = TRUE)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(quantile(x, probs = seq(0, 1, 0.1)))
	summary <- c((i-1), summary, mean(x))
	outtable[i, 1:length(summary)]<-summary
}

write.csv(outtable,file=paste0(basedir,"/percentileE28.csv"), row.names = FALSE)