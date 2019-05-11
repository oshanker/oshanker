#basedir<-"../oldriemann/out/gzetaE28"
basedir<-"../oldriemann/out/gzetaE12"
#outfile<-"/percentileE28.csv"
outfile<-"/percentile_calc.csv"
#infile<-"/gzeta6.csv"
infile<-"/gzeta_calc6.csv"


#lnfactor<-sqrt(log(1.0E12/(2*pi)))
df <- read.csv(paste0(basedir,infile), header = FALSE);
outtable<-matrix(nrow =  length(df), ncol = 15, byrow = TRUE)
quantileValues = c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1.0)
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(quantile(x, probs = quantileValues))
	summary <- c((i-1), summary, mean(x))
	outtable[i, 1:length(summary)]<-summary
}

write.csv(outtable,file=paste0(basedir, outfile), row.names = FALSE)
