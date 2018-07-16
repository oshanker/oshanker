library(splines)
 x=1:10
 y=x^2+1
df=data.frame(x,y)

print(summary(fm1 <- lm(y ~ bs(x, knots=c(2,3,5)), data = df)))
ht <- seq(3, 6, length.out = 4)
test=predict (fm1 ,newdata =list(x =ht),se=T)
print(ht)
print(test)