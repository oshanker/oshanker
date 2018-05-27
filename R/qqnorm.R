df <- read.csv("../oldriemann/data/zetaE12.csv");
ls(df)
x<-as.matrix(df[2]);
qqnorm(x, main = "Normal Q-Q Plot for Z(t) at Gram Points",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles"); qqline(x, col = 2)