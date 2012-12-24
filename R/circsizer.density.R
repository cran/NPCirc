circsizer.density <- function(x,NU,ngrid=250,alpha=0.05,B=500,type=3,raw.data=FALSE,log.scale=TRUE,
zero=pi/2,clockwise=TRUE,title=NULL,labels=NULL,label.pos=NULL,rad.pos=NULL){

	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	x <- x[!is.na(x)]
	n <- length(x)
	if (sum(is.na(x))>0) warning("Missing values were removed")
	if (n==0) stop("No observations (at least after removing missing values)")
	if (any(x<0) | any(x>2*pi)) stop("The sample of angles 'x' must be in radians between 0 and 2*pi")
	if (!is.numeric(NU)) stop("Argument 'NU' must be numeric")
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

	if (type==5){
		label.pos <- seq(pi/12,23*pi/12,by=pi/6)
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
	
	rad<-seq(0,2*pi,length=250)
	nu <- 1
	if ((!log.scale & NU[1]==0) | (log.scale & NU[nNU]==0)) { 
		radial.plot(c(rep(0,250),rep(stNU,250)),c(rad,rad),rp.type="p",poly.col=6,line.col=6,radial.lim=radlim,add=T,start=zero,clockwise=clockwise)
		nu<-2
	}
	t <- seq(0,2*pi,length=ngrid)
	t_x <- outer(t,x,"-")
	
	for (i in nu:nNU){

		exptx <- exp(cos(t_x)*NU[i])
		exptxsin <- exptx*sin(t_x)  
		pibessel <- 2*pi*besselI(NU[i],0)

		# Derivative estimate
		deriv_fhat <- -NU[i]*rowMeans(exptxsin)/pibessel 

		# Derivative of the von Mises kernel
		deriv_kernel <- -NU[i]*exptxsin/pibessel

		if (sum(is.na(deriv_fhat))>0 | sum(is.na(deriv_kernel))>0 | sum(deriv_kernel==Inf)>0){
			stop("Values of the smoothing parameter too large")
		}

		# Standard deviation	
		dt <- sqrt((rowMeans(deriv_kernel^2)-rowMeans(deriv_kernel)^2)/n)

		# Effective Sample Size
		ESS <-rowSums(exptx)/exp(NU[i]) 
		
		# Bootstrap-t confidence interval	
		zstar<-matrix(NA,B,ngrid)
		for (b in 1:B){
			xstar<- sample(x,replace=T)
			t_xstar <- outer(t,xstar,"-")
			exptxstar <- exp(cos(t_xstar)*NU[i])
			exptxstarsin <- exptxstar*sin(t_xstar)
			deriv_fhat_xstar <- -NU[i]*rowMeans(exptxstarsin)/pibessel
			deriv_kernel <- -NU[i]*exptxstarsin/pibessel 
			dtest <- sqrt((rowMeans(deriv_kernel^2)-rowMeans(deriv_kernel)^2)/n)
			zstar[b,]<-(deriv_fhat_xstar-deriv_fhat)/dtest
		}
		quant <- apply(zstar, 2, function(z) quantile(z, probs=c(alpha/2,1-alpha/2),type=1))
		li<-deriv_fhat-quant[2,]*dt
		ls<-deriv_fhat-quant[1,]*dt

		# Add the color rings		
		color<-rep(6,ngrid)
		color[which(li>0)]<-"darkblue"
		color[which(ls<0)]<-"red"
		color[which(ESS<5)]<-"grey"
		if (log.scale) NUi<-NUplot[i]-stNU
		else NUi<-ifelse(NUplot[i]-stNU<0,0,NUplot[i]-stNU)
		NUs<-NUplot[i]+stNU
		for (j in 1:ngrid){
			stt<-(t[2]-t[1])/2
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
	# Add the raw data
	if (raw.data){
		bins <- n
		x <- -x + pi/2
		x <- x%%(2 * pi)
		x[x >= 2 * pi] <- 2 * pi - 4 * .Machine$double.eps
		arc <- (2 * pi)/bins
		pos.bins <- (1 - 1/2) * arc - arc/2
		breaks <- seq(0, 2 * pi, length.out = (bins + 1))
		bins.count <- hist.default(x, breaks = breaks, plot = FALSE, right = TRUE)$counts
		mids <- seq(arc/2, 2 * pi - pi/bins, length = bins) + pos.bins[1]
		for (i in 1:bins) {
			if (bins.count[i] != 0) {
				for (j in 0:(bins.count[i] - 1)) {
					r <- 1 + j * 0.025
                  		z <- r * cos(mids[i])+1
                  		y <- r * sin(mids[i])
					points.default(z*maxlength-maxlength,y*maxlength, pch = 16, cex = 1)
				}
			}
		}
	}
}
