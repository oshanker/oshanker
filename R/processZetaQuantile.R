basedir<-"../oldriemann/out/gzetaE12"
intable <- read.csv(paste0(basedir,"/percentile_calc.csv"),header=TRUE);
meanValues<- c(1:13)
slopes<- c(1:13)
for (i in 2:13) {
	x<-as.matrix(intable[i]);
	x<-as.vector(x);
	print(x)
	print(paste('is',i))
	meanValues[i]=mean(x) 
	slopes[i]=(x[3]-x[5])/2
}
    conv = format( meanValues, scientific = FALSE, drop0trailing = TRUE, digits = 4, width =8)
#    conv = as.numeric(conv)
    print(conv, quote = FALSE)
#print(meanValues)
print(slopes)