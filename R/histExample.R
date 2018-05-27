df <- read.csv("../oldriemann/data/zetaE12.csv");
ls(df)
x<-as.matrix(df[2]);
x<-as.vector(x);
print("Data for ../oldriemann/data/zetaE12.csv")
print(ls.str())
print("Stats for ../oldriemann/data/zetaE12.csv")
print(quantile(x))
print(summary(x))
y<-1:length(x)
even<-y[y%%2==0]
print("Stats for even")
print(summary(x[even]))
odd<-y[y%%2==1]
print("Stats for odd")
print(summary(x[odd]))
cross=x[odd]*x[even]
print("Stats for cross-product")
print(summary(cross))
b=seq(-100,25,5)
b=c(-4000,b,12000)
h=hist(cross,breaks=b,xlim=c(-100,25),freq=TRUE)

