basedir<-"../oldriemann/data/gzetaE12"
infile<-"/quantileFitted_calc6.csv"
basedir1<-"../oldriemann/data/gzetaE28"
infile1<-"/quantileFittedE28.csv"
ranges=3:13
plotindex<-2
intable <- read.csv(paste0(basedir,infile),header=FALSE);
	x<-as.matrix(intable[1,ranges]);
	x<-as.vector(x);

	y0<-as.matrix(intable[plotindex,ranges]);
	y0<-as.vector(y0);
	
y<-matrix(nrow=length(x),ncol=4)
	print(dim(y))
	y[1:length(x),1]<-y0

intable1 <- read.csv(paste0(basedir1,infile1),header=FALSE);
	y1<-as.matrix(intable1[plotindex,ranges]);
	y1<-as.vector(y1);
	
	y[1:length(x),2]<-y1

print(x)
print(y)

f=0.1
den=0.25-f^2
a=(y[7,1:2]-y[5,1:2])*den/(2*f)

f=x[11]-0.5
den=0.25-f^2
b=(y[11,1:2]*den-a*f)/(f*f*f)

print(paste(c("a12","a28"),a))
print(paste(c("b12","b28"),b))

for (i in 1:length(x)) {
   xx=x[i]-0.5
   den=0.25-xx^2
   y[i,3:4]=(a*xx+b*xx*xx*xx)/den
}

matplot(x, y, type = "b", xaxt = "n",#
        main = "Quantile",#
        pch = "1*ab", 
        ylab = expression("mean Z"), # only 1st is taken#
        xlab = expression("f")#
        )#
axis(1)
text(0.84,3.9, "E28")
text(0.9,2.8, "E12")
grid()
