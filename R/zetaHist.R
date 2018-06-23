df <- read.csv("../oldriemann/out/gzetaE12/gzeta.csv");
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	print(i-1)
	summary <- summary(x)
    print( summary)
}
