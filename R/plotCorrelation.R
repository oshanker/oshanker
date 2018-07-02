df <- read.csv("../oldriemann/out/gzetaCorrelation/gzeta8.csv", header = FALSE);
ls(df)
correlation = as.vector(df[9])
xx = seq(1,16) 
matplot(xx,correlation,type = "l")