modal.clustering<-function(x,bw=NULL,tol=0.0000000000001, 
                           plot=TRUE,shrink = 1, pch=16, col.points="grey",
                           stack=TRUE, labels=NULL, control.circular = list()){
  
  
  
  data <- x
  name <- deparse(substitute(x))
  
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric")
  if (is.circular(x)){
    datacircularp <- circularp(x)
  }else {
    datacircularp <- list(type = "angles", units = "radians", 
                          template = "none", modulo = "asis", zero = 0, rotation = "counter")
  }
  dc <- control.circular
  if (is.null(dc$type)) 
    dc$type <- datacircularp$type
  if (is.null(dc$units)) 
    dc$units <- datacircularp$units
  if (is.null(dc$template)) 
    dc$template <- datacircularp$template
  if (is.null(dc$modulo)) 
    dc$modulo <- datacircularp$modulo
  if (is.null(dc$zero)) 
    dc$zero <- datacircularp$zero
  if (is.null(dc$rotation)) 
    dc$rotation <- datacircularp$rotation
  if (dc$template=="clock12" | dc$template=="clock24"){
    stop("function not yet supported for template='clock12' or template='clock24")
  }
  if (dc$units=="hours"){stop("function not yet supported for units='hours'")}
  
  
  if (dc$modulo == "pi"){x<-2*x}
  x <- conversion.circular(x, units = "radians", zero = 0, 
                           rotation = "counter",modulo="asis")
  #attr(x, "class") <- attr(x, "circularp") <- NULL
  
  #x <- as.vector(x)
  x.na <- is.na(x)
  if (any(x.na)) {
    x <- x[!x.na]
    warning("x contains missing values")
  }
  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
  }
  
  
  if (!is.null(bw)&!is.numeric(bw)) {stop("argument 'bw' must be numeric")}
  if (!is.null(bw)){if(!is.finite(bw)) {stop("non-finite `bw'")}}
  if (!is.null(bw)){if(bw < 0) {stop("`bw' is not positive")}}
  
  if (!is.null(tol)&!is.numeric(tol)) {stop("argument 'tol' must be numeric")}
  if (!is.null(tol)){if(!is.finite(tol)) {stop("non-finite `tol'")}}
  if (!is.null(tol)){if(tol <= 0) {stop("`tol' is not positive")}}
  
  pch <- check_pch(pch, "pch")
  
  if(!is.logical(plot)){stop("argument 'plot' must be logical")}
  if(!is.logical(stack)){stop("argument 'stack' must be logical")}
  
  if (length(col.points)!= 1){stop("'col.points' must be a single value (number or color name)")}
  if (!is.numeric(col.points) & !is.character(col.points)){stop("'col.points' must be either a numeric index or a color name (character)")}
  
  if (!is.numeric(shrink) || length(shrink) != 1 || shrink <= 0) {
    stop("'shrink' must be a single positive number")
  }
  
  if(!is.null(labels)){if(!is.character(labels)){stop("'labels' must be character" )}}
  
  
  n<-length(x)
  estimation<-mode_estimation(x,bw=bw,tol=tol)
  modes<-estimation$modes
  bw<-estimation$bw
  clust<-estimation$clust
  cols<-estimation$cols
  
  modes<-circular(modes,type="angles",units="radians",modulo="2pi",zero=0,rotation="counter")
  
  modes_transformed<-conversion.circular(modes,type=dc$type,
                                         units=dc$units,template=dc$template,
                                         zero=dc$zero,
                                         rotation=dc$rotation,
                                         modulo="2pi")
  
  
  if (dc$modulo == "pi"){
    modes_transformed<-modes_transformed/2;
    #modes_transformed<-conversion.circular(modes_transformed,modulo=dc$modulo)
  }
  
  if(plot==TRUE){
    
    if(dc$modulo=="pi"){
      if(dc$units=="radians"){
        opo<-pi
      }else if(dc$units=="degrees"){
        opo<-180
      }
      
    }
    
    plot(data,col=col.points,shrink=shrink,pch=pch,
         axes=FALSE,stack=stack)
    if(dc$modulo=="pi"){
      points.circular(data+opo,stack=TRUE,col=col.points)
    }
    # add data clustered by color
    for(k in 1:n){
      ll<-data[k]
      segments.circular(ll,col=cols[k],lwd=2)
    }
    if(dc$modulo=="pi"){
      for(k in 1:n){
        ll<-data[k]+opo
        segments.circular(ll,col=cols[k],lwd=2)
      }
    }
    
    # add modes
    for (j in 1:length(modes_transformed)){
      cc<-modes_transformed[j]
      segments.circular(cc,col=1,lwd=2)
    }
    if(dc$modulo=="pi"){
      for (j in 1:length(modes_transformed)){
        cc<-modes_transformed[j]+opo
        segments.circular(cc,col=1,lwd=2)
      }
    }
    if(is.null(labels)){
      axis.circular(units=dc$units,template=dc$template,zero=dc$zero,
                    rotation=dc$rotation,tcl.text = -0.15)
    }else{
      k<-length(labels)
      axis.circular(seq(0, (2*k - 2)*pi/k, length.out = k),labels,tcl.text = -0.15)
    }
    
  }
  return(list(modes=modes_transformed,bw=bw,clust=clust))
}
