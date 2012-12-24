circsizer.regression <- function(x,y,NU,ngrid=150,alpha=0.05,B=500,B2=250,type=3,log.scale=TRUE,
zero=pi/2,clockwise=TRUE,title=NULL,labels=NULL,label.pos=NULL,rad.pos=NULL){

	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.numeric(y)) stop("argument 'y' must be numeric")
	if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
	nax <- is.na(x)
	nay <- is.na(y)
	x<-x[!nax & !nay]
	y<-y[!nax & !nay]
	if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.")
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (any(x<0) | any(x>2*pi)) stop("The sample of angles 'x' must be in radians between 0 and 2*pi")
	if (!is.numeric(NU)) stop("argument 'NU' must be numeric")
	if (length(NU)==1) stop("Length of argument 'NU' must be at least 2")
	if (any(NU<0)) warning("Values of vector 'NU' must be positive. Values smaller than 0 were removed")
	NU <- sort(NU[NU>=0]) 
	if (!log.scale){ 
		NUeq <- seq(min(NU),max(NU),length=length(NU))
		if (all(NU==NUeq)==FALSE){
			NU <- NUeq
			warning ("Values of 'NU' are not equally spaced. CircSiZer is computed for a equally spaced grid between 
			the minimum and maximum values of 'NU'")
		}
		NUplot <- NU

	}else{ 
		NUplot <- seq(-log10(max(NU[NU!=0])),-log10(min(NU[NU!=0])),length=length(NU))
		NU <- 10^(-NUplot)
		if(length(unique(NU))==1){
			warning("Argument 'NU' is not valid. A default grid was used")
			NUplot<-seq(-2,1,by=0.1)
			NU<-10^(-NUplot)
		}
	}
	if (!is.numeric(ngrid)) stop("argument 'ngrid' must be numeric")
	if (ngrid<=0){
		warning("'ngrid' must be positive integer. Default value of 'ngrid' was used")
		ngrid=250
	}
	if (!is.numeric(alpha)) stop("argument 'alpha' must be numeric")
	if (alpha<0 | alpha>1){
		warning("'alpha' must be in the interval [0,1]. Default value of 'alpha' was used")
		alpha=0.05
	}
	if (!is.numeric(B)) stop("argument 'B' must be numeric")
	if (B<=0){
		warning("'B' must be positive. Default value of 'B' was used")
		B=500
	}
	if (!is.numeric(B2)) stop("argument 'B' must be numeric")
	if (B2<=0){
		warning("'B2' must be positive. Default value of 'B2' was used")
		B2=250
	}
	if (is.null(labels)){
		if (!any(type==1:5)){ 
			type=3
			warning("Value specified for argument 'type' is not valid. 'type=3' was used")
		}
		label.pos <- rad.pos <- seq(0,7*pi/4,by=pi/4)
	}else{ 
		if(length(labels)!=length(label.pos)){
			warning("Arguments 'labels' and 'label.pos' are not the same length. 'type=3' was used")
			type=3
			label.pos <- rad.pos <- seq(0,7*pi/4,by=pi/4)
		}else{
			type=0
			if (is.null(rad.pos)) rad.pos <- seq(0,7*pi/4,by=pi/4)
		}
	}

	if (type==5){label.pos <- seq(pi/12,23*pi/12,by=pi/6)
			 rad.pos <- seq(0,11*pi/6,by=pi/6) 
	}
	if (type==1) labels <- c("N","NE","E","SE","S","SW","W","NW")
	if (type==2) labels <- c("0h","3h","6h","9h","12h","15h","18h","21h")
	if (type==3) labels <- c(expression(0),expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi),
		                   expression(5*pi/4),expression(3*pi/2),expression(7*pi/4))
	if (type==4) labels <- c("0","45","90","135","180","225","270","315")
	if (type==5) labels <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Ago","Sep","Oct","Nov","Dec")

	# Create the plot
	nNU <- length(NUplot) 
	stNU <- (NUplot[2]-NUplot[1])/2 
	radlim<- c((NUplot[1]-stNU),NUplot[nNU]+stNU)
	radlim[1]<-radlim[1]-diff(radlim)*0.2
	grid.pos<-pretty(radlim)
	st<-grid.pos[2]-grid.pos[1]
	radial.plot(0,rp.type="s",labels="",point.col="white",radial.lim=radlim, show.grid=F,show.radial.grid=F,
     	radial.labels="",show.grid.labels=F,start=zero,clockwise=clockwise,main=title)   
	# Add the labels
	label.prop=1.12
	maxlength=diff(radlim)
	if (clockwise) label.pos <- -label.pos + zero
	else label.pos <- label.pos + zero
	xpos <- cos(label.pos) * maxlength * label.prop
	ypos <- sin(label.pos) * maxlength * label.prop
	text(xpos, ypos, labels, cex = par("cex.axis"))
	# Plot the arrow
	arrow.pos <- seq(pi/4-0.13, pi/4+0.13, length=250)  
	xarrowpos <- cos(arrow.pos) * maxlength * 1.28
	yarrowpos <- sin(arrow.pos) * maxlength * 1.28
	points.default(xarrowpos,yarrowpos,type="l",lwd=2)
	if (clockwise) Arrows(xarrowpos[1], yarrowpos[1], xarrowpos[1]+0.015, yarrowpos[1]-0.015, type="simple", arr.length=0.5, segment=FALSE)
	else Arrows(xarrowpos[250], yarrowpos[250], xarrowpos[250]-0.015, yarrowpos[250]+0.015, type="simple", arr.length=0.5, segment=FALSE)
	
	# Auxiliary function that computes the estimator of the derivate of the regression function
	der.cl<-function(x,y,nu){
		x_t <- matrix(x,nrow=n)%*%v1t - v1xt
		m<-sin(x_t)
		kxt<-exp(nu*cos(x_t))
		sn1<-colSums(kxt*m)   
		sn2<-colSums(kxt*m^2) 
		sn3<-colSums(kxt)     
		titawy1<-colSums(kxt*y)   
		titawy2<-colSums(kxt*m*y)
		beta1<- -sn1*titawy1+sn3*titawy2
		det<-sn3*sn2-sn1^2 
		return(beta1/det)
	}

	t <- seq(0,2*pi,length=ngrid)
	v1t<-rep(1,ngrid)
	v1xt<-matrix(rep(1,n),nrow=n)%*%t

	nu <- 1
	rad<-seq(0,2*pi,length=250)
	if ((!log.scale & NU[1]==0) | (log.scale & NU[nNU]==0)) { 
		radial.plot(c(rep(0,250),rep(stNU,250)),c(rad,rad),rp.type="p",poly.col=6,line.col=6,radial.lim=radlim,add=T,start=zero,clockwise=clockwise)
		nu<-2
	}
	for (i in nu:nNU){

		# Estimate of the derivate of the regression function
		b <- der.cl(x,y,NU[i]) 
		if (sum(is.na(b))>0){
			stop("Values of the smoothing parameter too large")
		}
		
		# Standard deviation
		Xstar<-Ystar<-matrix(0,B,n)
		der<-matrix(0,B,ngrid)
		ind<-matrix(sample(n,n*B,replace=T),n,B)
		for (k in 1:B){
			indk<-ind[,k]
			xstar<-x[indk]
			ystar<-y[indk]
			der[k,] <- der.cl(xstar,ystar,NU[i]) 
			Xstar[k,]<-xstar
			Ystar[k,]<-ystar
		}

		dt<-sqrt(colMeans(der^2)-colMeans(der)^2)
		
		# Effective Sample Size
		exptx<-exp(NU[i]*cos(outer(t,x,"-")))
		ESS <-rowSums(exptx)/exp(NU[i])
		
		# Bootstrap-t confidence interval		
		zstar<-matrix(0,B,ngrid)
		derstarstar<-matrix(0,B2,ngrid)
		for (l in 1:B){
			ind<-matrix(sample(n,n*B2,replace=T),n,B2)
			for (k in 1:B2){
				indk<-ind[,k]
				xstarstar<-Xstar[l,][indk]
				ystarstar<-Ystar[l,][indk]
				derstarstar[k,] <- der.cl(xstarstar,ystarstar,NU[i])
			}
			dtest<-sqrt(colMeans(derstarstar^2)-colMeans(derstarstar)^2)
			zstar[l,]<-(der[l,]-b)/dtest
		}
		quant <- apply(zstar, 2, function(z) quantile(z, probs=c(alpha/2,1-alpha/2),type=1))
		li<-b-quant[2,]*dt
		ls<-b-quant[1,]*dt
		
		# Add the color rings
		color<-rep(6,ngrid)
		color[which(li>0)]<-"darkblue"
		color[which(ls<0)]<-"red"
		color[which(ESS<5)]<-"grey"
		for (j in 1:ngrid){
			stt<-(t[2]-t[1])/2
			if (log.scale) NUi<-NUplot[i]-stNU
			else NUi<-ifelse(NUplot[i]-stNU<0,0,NUplot[i]-stNU)
			NUs<-NUplot[i]+stNU
			ti<-t[j]-stt
			ts<-t[j]+stt
			v1<-c(ti,NUi)
			v2<-c(ti,NUs)
			v3<-c(ts,NUi)
			v4<-c(ts,NUs)
			v1v3<-seq(v1[1],v3[1],by=0.01)
			v4v2<-sort(v1v3,decreasing=T)
			radial.plot(c(v1[2],rep(NUi,length(v1v3)),v3[2],v4[2],rep(NUs,length(v4v2)),v2[2]),c(v1[1],v1v3,v3[1],v4[1],v4v2,v2[1]),rp.type="p",poly.col=color[j],line.col=color[j],radial.lim=radlim,add=T,start=zero,clockwise=clockwise)
		}
	}
	# Add the rings
	xcir <- grid.pos
	ncir <- length(xcir)
	ind<-c(which(xcir>NUplot[nNU]+stNU),which(xcir<NUplot[1]-stNU))
	if(length(ind)>0){xcir<-xcir[-ind]; ncir<-ncir-length(ind)}
	for(i in 1:ncir){radial.plot(rep(xcir[i],250),rad,rp.type="p",radial.lim=radlim,add=T)} 
	# Add the radius
	segments(rep(0,9),rep(0,9),cos(rad.pos)*maxlength,sin(rad.pos)*maxlength) 
	# Add the values of the smoothing parameter along the radius
	xpos <- xcir - radlim[1]
	ypos <- rep(0, ncir)
	boxed.labels(xpos,ypos,as.character(xcir),border = FALSE,cex=0.7) 
}