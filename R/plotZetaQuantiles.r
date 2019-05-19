basedir<-c("../oldriemann/data/gzetaE12","../oldriemann/data/gzetaE28")
infile<-c("/quantileFitted_calc6.csv", "/quantileFittedE28.csv")
#basedir1<-"../oldriemann/data/gzetaE28"
#infile1<-"/quantileFittedE28.csv"
ranges=3:13
plotindex<-2
y<-matrix(nrow=length(ranges),ncol=4)
	print(dim(y))

intable <- read.csv(paste0(basedir[1],infile[1]),header=FALSE);
	f<-as.matrix(intable[1,ranges]);
	f<-as.vector(x);

	y0<-as.matrix(intable[plotindex,ranges]);
	y0<-as.vector(y0);
	
	y[1:length(x),1]<-y0

intable1 <- read.csv(paste0(basedir[2],infile[2]),header=FALSE);
	y1<-as.matrix(intable1[plotindex,ranges]);
	y1<-as.vector(y1);
	
	y[1:length(x),2]<-y1

print(f)
print(y)

ind=c(rep(0,10),rep(1,10))
xxvec=f-0.5
xx2=xxvec ^2
	xx2trun=c(xx2[1:5],xx2[7:11],xx2[1:5],xx2[7:11])
	xxtrun=c(xxvec[1:5],xxvec[7:11],xxvec[1:5],xxvec[7:11])
	q=c(y[1:5,1],y[7:11,1],y[1:5,2],y[7:11,2])
	xxq=xxtrun/q
	lm.fit2=lm(xx2trun~ ind+xxq )
	store=coef(lm.fit2)
	print(summary (lm.fit2))
	print(store)


poles=c(0.263, 0.250)
#poles=c(store[1],store[1]+store[2])
for (i in 1:2) {
	pp=rep(poles[i],11)
    denvec=pp-xx2
	yy0=y[1:length(xxvec),i]*denvec
	lm.fit2=lm(yy0~ xxvec )
	#lm.fit2=lm(yy0~ xxvec +I(xxvec ^3))
	print(summary (lm.fit2))
	fit0=coef(lm.fit2)
	y[1:length(xxvec),i+2]=(fit0[2]*xxvec)/denvec
}



matplot(f, y, type = "b", xaxt = "n",#
        main = "Quantile",#
        pch = "1*ab", 
        ylab = expression("mean Z"), # only 1st is taken#
        xlab = expression("f")#
        )#
axis(1)
text(0.84,3.9, "E28")
text(0.9,2.8, "E12")
legend("topleft", c("E12","E28","fit E12","fit E28"),pch = "1*ab")
grid()
