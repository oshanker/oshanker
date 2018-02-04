df <- read.csv("../oldriemann/data/zetaE12.csv");
ls(df)
x<-as.matrix(df[2]);
x<-as.vector(x);
print(ls.str())
print(quantile(x))
print(summary(x))
y<-1:length(x)
even<-y[y%%2==0]
print(summary(x[even]))
odd<-y[y%%2==1]
print(summary(x[odd]))