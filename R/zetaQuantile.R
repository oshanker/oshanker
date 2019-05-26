basedir<-c("../oldriemann/out/gzetaE12", "../oldriemann/out/gzetaE28")
outfile<-c("/percentile_calc.csv","/percentileE28.csv")
infile<-c("/gzeta_calc6.csv","/gzeta6.csv")
runFor=2

#lnfactor<-sqrt(log(1.0E12/(2*pi)))
df <- read.csv(paste0(basedir[runFor],infile[runFor]), header = FALSE);
quantileValues = c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1.0)
outtable<-matrix(nrow =  length(df), ncol = 3+length(quantileValues), byrow = TRUE)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(quantile(x, probs = quantileValues))
	summary <- c((i-1), summary, mean(x),sd(x))
	outtable[i, 1:length(summary)]<-summary
}

write.csv(outtable,file=paste0(basedir[runFor], outfile[runFor]), row.names = FALSE)
