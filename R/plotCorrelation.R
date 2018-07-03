df <- read.csv("../oldriemann/out/gzetaCorrelation/gzeta8.csv", header = FALSE);
ls(df)
i = c(1,9)
correlation = as.vector(df[i])
xx = seq(1,16) 
matplot(xx,correlation,type = "b", tck=1)