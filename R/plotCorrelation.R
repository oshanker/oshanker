df <- read.csv("../oldriemann/out/gzetaCorrelation/gzeta8.csv", header = FALSE);
ls(df)
i = c(1,9)
transposed= t(as.matrix(df))
correlation = (transposed[,i])
xx = t(seq(1,16)) 
xx = c(xx, xx)
dim(xx) = c(16,2)
matplot(xx,correlation,type = "b", tck=1)