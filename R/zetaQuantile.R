basedir<-c("../oldriemann/out/gzetaE12", "../oldriemann/out/gzetaE28")
outfile<-c("/percentile_calc.csv","/percentileE28.csv")
sdoutfile<-c("/sd_calc.csv","/sdE28.csv")
infile<-c("/gzeta_calc6.csv","/gzeta6.csv")
runFor=1

#lnfactor<-sqrt(log(1.0E12/(2*pi)))
df <- read.csv(paste0(basedir[runFor],infile[runFor]), header = FALSE);
quantileValues = c(0,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,1.0)
outtable<-matrix(nrow =  1+length(df), ncol = 1+length(quantileValues), byrow = TRUE)
sdouttable<-matrix(nrow =  3, ncol = length(df), byrow = TRUE)
outtable[1, 1:dim(outtable)[2]]<-c(0,quantileValues)
	
for (iangle in 1:length(df)) {
	x<-as.matrix(df[iangle]);
	x<-as.vector(x);
	summary <- as.vector(quantile(x, probs = quantileValues))
	summary <- c((iangle-1), summary)
	sdsummary <- c((iangle-1),  mean(x),sd(x))
	outtable[iangle+1, 1:length(summary)]<-summary
	sdouttable[1:length(sdsummary),iangle]<-sdsummary
}

write.csv(outtable,file=paste0(basedir[runFor], outfile[runFor]), row.names = FALSE)
write.csv(sdouttable,file=paste0(basedir[runFor], sdoutfile[runFor]), row.names = FALSE)
