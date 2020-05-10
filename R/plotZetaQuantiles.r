basedir<-c("../oldriemann/data/gzetaE12","../oldriemann/data/gzetaE28")
infile<-c("/quantileFitted_calc6.csv", "/quantileFittedE28.csv")
ranges=3:13
# cover f from 0.1 to 0.9
plotindex <- 2
# 3 : slope; 2 : intercept
y<-matrix(nrow=length(ranges),ncol=4)
print(dim(y))

intable <- read.csv(paste0(basedir[1],infile[1]),header=FALSE);
f<-as.matrix(intable[1,ranges]);
f<-as.vector(f);

assignDependent <- function(intable, plotindex, ranges) {
  y0<-as.matrix(intable[plotindex,ranges]);
  y0<-as.vector(y0);
  
  return(y0)
}

intable1 <- read.csv(paste0(basedir[2],infile[2]),header=FALSE);

xxvec=f-0.5
xx2=xxvec ^2
runpolefit=FALSE
if(runpolefit) {
  ind=c(rep(0,10),rep(1,10))
	xx2trun=c(xx2[1:5],xx2[7:11],xx2[1:5],xx2[7:11])
	xxtrun=c(xxvec[1:5],xxvec[7:11],xxvec[1:5],xxvec[7:11])
	q=c(y[1:5,1],y[7:11,1],y[1:5,2],y[7:11,2])
	xxq=xxtrun/q
	lm.fit2=lm(xx2trun~ ind+xxq )
	store=coef(lm.fit2)
	print(summary (lm.fit2))
	print(store)
}

legendposition="topleft"

print(f)
for (plotindex in range(2,3)) {
  
  y[1:length(f),1] <- assignDependent(intable, plotindex, ranges)
  y[1:length(f),2] <- assignDependent(intable1, plotindex, ranges)
  print(y)
  
  if(plotindex==2){
  	a=1:2
  	poles2=c(0.29, 0.252)
  	#poles=c(store[1],store[1]+store[2])
  	for (i in 1:2) {
  		pp=rep(poles2[i],length(xxvec))
  	    denvec=pp-xx2
  		yy0=y[1:length(xxvec),i]*denvec
  		lm.fit2=lm(yy0~ xxvec )
  		#lm.fit2=lm(yy0~ xxvec +I(xxvec ^3))
  		print(summary (lm.fit2))
  		fit0=coef(lm.fit2)
  		a[i]=(fit0[2])
  		y[1:length(xxvec),i+2]=(fit0[2]*xxvec)/denvec
  	}
  	rootcheck=(a[1]*poles2[2]-a[2]*poles2[1])/(a[1]-a[2])
  	print(paste('rootcheck', sqrt(rootcheck) ) )
  	label = 'intercept'
  } else if(plotindex==3){
  	legendposition="center"
  	poles2=c(0.41, 0.45)
  	for (i in 1:2) {
  		pp=rep(poles2[i],length(xxvec))
  	    denvec=pp-xx2
  		yy0=y[1:length(xxvec),i]*denvec
  		lm.fit2=lm(yy0~I(xx2))
  		print(summary (lm.fit2))
  		fit0=coef(lm.fit2)
  		print(fit0)
  		y[1:length(xxvec),i+2]=(fit0[1] + fit0[2]*xx2)/denvec
  		label = 'slope'
  	}
  }
  
  symbols="12a*"
  matplot(f, y, type = "b", xaxt = "n",#
          main = "Quantile",#
          pch = symbols, 
          ylab = label, # only 1st is taken#
          xlab = expression("f")#
          )#
  axis(1)
  text(0.84, 3.9, "E28")
  text(0.9, 2.8, "E12")
  legend(legendposition, c("E12","E28","fit E12","fit E28"),pch = symbols)
  grid()
  opt = readline(paste("completed plotindex", plotindex ))
  if(opt == "q"){
    breakq
  }
}