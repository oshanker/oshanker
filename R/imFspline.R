library(splines)

df <- read.csv("../oldriemann/out/imFmidGramE12.csv", header = FALSE);
x=as.vector(df[1])
breaks=seq(3, length(x[,1])-2, length.out = length(x[,1]) - 4)
fm1 <- lm(V2 ~ bs(V1, knots=df[breaks,1]), data = df)
ht <- df[1:length(x[,1])-1,1]+0.5
test=predict (fm1 ,newdata =list(V1 =ht),se=F)

figure = function() {
#plot(ht ,test ,col =" gray ", type="l")
#lines(df$V1 ,df$V2)
xx=c(df$V1 ,ht)
yy=c(df$V2, test)
plot(xx,yy, type="p", tck=1, pch=".")
lines(ht,test, col="green")
lines(df$V1 ,df$V2, col="red")
}

#figure()
out = data.frame(ht+0.5, test)
write.csv(out, row.names = F, file = "../oldriemann/out/imFGramE12.csv")