rm(list=ls())
myfunc <- function(i, df){
#for (i in 1:len1) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(summary(x))
#	print(summary(x))
	summary <- c((i), summary)
    conv = format( summary, scientific = FALSE, drop0trailing = TRUE, digits = 4, width =8)
#    conv = as.numeric(conv)
    print(conv, quote = FALSE)
#}


lower = -1
upper = 1       
b=seq(lower,upper,0.1)
yy = range(x)
b=c(yy[1],b,yy[2])
h=hist(x,breaks=b,xlim=c(lower,upper),xlab='Z',
   main = paste("Histogram of Z for " , (i-1),  expression(pi), "/12"))
}

df <- read.csv("../oldriemann/out/gzetaE12/gzeta_calc12.csv", header = FALSE);
len1=length(df)
myfunc(13,df)
scan("")
myfunc(1,df)