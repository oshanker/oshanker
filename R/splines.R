library(splines)
 x=1:10
 y=c(2.05,   5.1,  10.2,  17.0,  26.1,  37.3,  50.2,  65.4,  82.7, 101.5)
df=data.frame(x,y)

print(summary(fm1 <- lm(y ~ bs(x, knots=3:8), data = df)))
ht <- seq(3, 6, length.out = 4)
test=predict (fm1 ,newdata =list(x =ht),se=F)
print(ht)
print(test)