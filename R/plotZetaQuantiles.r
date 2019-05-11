basedir<-"../oldriemann/data/gzetaE12"
infile<-"/quantileFitted_calc6.csv"
basedir1<-"../oldriemann/out/gzetaE28"
infile1<-"/quantileFittedE28.csv"
intable <- read.csv(paste0(basedir,infile),header=FALSE);
	x<-as.matrix(intable[1,3:13]);
	x<-as.vector(x);

	y0<-as.matrix(intable[2,3:13]);
	y0<-as.vector(y0);
	
y<-matrix(nrow=length(x),ncol=2)
	print(dim(y))
	y[1:length(x),1]<-y0

intable1 <- read.csv(paste0(basedir1,infile1),header=FALSE);
	y1<-as.matrix(intable1[2,3:13]);
	y1<-as.vector(y1);
	
	y[1:length(x),2]<-y1

print(x)
print(y)

matplot(x, y, type = "l", xaxt = "n",#
        main = "Quantile",#
        ylab = expression("mean Z"), # only 1st is taken#
        xlab = expression("Quantile ")#
        )#
axis(1)
text(0.9,3.9, "E28")
text(0.9,2.8, "E12")
