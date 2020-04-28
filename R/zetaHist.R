#input for E28 from math.MoreGSeriesTest.testInterpolate()
df <- read.csv("../oldriemann/out/gzetaE28/gzeta6.csv", header = FALSE);
for (i in 1:length(df)) {
	x<-as.matrix(df[i]);
	x<-as.vector(x);
	summary <- as.vector(summary(x))
#	print(summary(x))
	summary <- c((i-1), summary)
    conv = format( summary, scientific = FALSE, drop0trailing = TRUE, digits = 4, width =8)
#    conv = as.numeric(conv)
    print(conv, quote = FALSE)
}
