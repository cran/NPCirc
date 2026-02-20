circ.modes<-function(x,bw.modes=NULL,bw.dens=NULL,plot=TRUE,col.modes=2,col.dens=1,
                     lwd.modes=2,lwd.dens=2,tol=0.0000000000001, 
                     labels=NULL, control.circular = list()){
  
  
  
  data <- x
  name <- deparse(substitute(x))
  
  if (!is.numeric(x)) 
    stop("argument 'x' must be numeric")
  if (is.circular(x)) {
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
  
  
  if (!is.null(bw.modes)&!is.numeric(bw.modes)) {stop("argument 'bw.modes' must be numeric")}
  if (!is.null(bw.modes)){if(!is.finite(bw.modes)) {stop("non-finite `bw.modes'")}}
  if (!is.null(bw.modes)){if(bw.modes < 0) {stop("`bw.modes' is not positive")}}
  
  if (!is.null(bw.dens)&!is.numeric(bw.dens)) {stop("argument 'bw.dens' must be numeric")}
  if (!is.null(bw.dens)){if(!is.finite(bw.dens)) {stop("non-finite `bw.dens'")}}
  if (!is.null(bw.dens)){if(bw.dens < 0) {stop("`bw.dens' is not positive")}}
  
  
  
  if(!is.logical(plot)){stop("argument 'plot' must be logical")}
  if (length(col.modes)!= 1){stop("'col.modes' must be a single value (number or color name)")}
  if (!is.numeric(col.modes) & !is.character(col.modes)){stop("'col.modes' must be either a numeric index or a color name (character)")}
  if (length(col.dens)!= 1){stop("'col.dens' must be a single value (number or color name)")}
  if (!is.numeric(col.dens) & !is.character(col.dens)){stop("'col.dens' must be either a numeric index or a color name (character)")}
  
  lwd.modes <- check_lwd(lwd.modes, "lwd.modes")
  lwd.dens  <- check_lwd(lwd.dens,  "lwd.dens")
  
  if (!is.null(tol)&!is.numeric(tol)) {stop("argument 'tol' must be numeric")}
  if (!is.null(tol)){if(!is.finite(tol)) {stop("non-finite `tol'")}}
  if (!is.null(tol)){if(tol <= 0) {stop("`tol' is not positive")}}
  
  if(!is.null(labels)){if(!is.character(labels)){stop("argument 'labels' must be character")}}
  
  estimation<-mode_estimation(x,bw=bw.modes,tol=tol)
  modes<-estimation$modes
  bw<-estimation$bw
  
  modes<-circular(modes,type="angles",units="radians",modulo="2pi",zero=0,rotation="counter")
  
  
  
  modes_transformed<-conversion.circular(modes,type=dc$type,
                                         units=dc$units,template=dc$template,
                                         zero=dc$zero,
                                         rotation=dc$rotation,
                                         modulo="2pi")
  
  
  if (dc$modulo == "pi"){modes_transformed<-modes_transformed/2}
  
  
  if(plot==TRUE){
    
    
    if(dc$modulo=="pi"){
      if(dc$units=="radians"){
        opo<-pi
      }else if(dc$units=="degrees"){
        opo<-180
      }
    }
    
    # kde for density
    if(is.null(bw.dens)){
      kappa<-bw.AA(x,deriv.order=0)
    }else{
      kappa<-bw.dens
    }
    dens_x<-density.circular(x,bw=kappa)
    dens_x_modes<-kern.den.circ(x,unique(round(modes,2)),bw=kappa)$y
    
    
    if(is.null(labels)){
      
      
      if(dc$template=="geographics"){
        label.pos=seq(0, (2*8 - 2)*pi/8, length.out = 8)
        labels=c("E","NE","N","NW","W","SW","S","SE")
      }else if(dc$template=="clock12"){
        label.pos=seq(0, (2*12 - 2)*pi/12, length.out = 12)
        labels=c(3,2,1,12:4)
      }else if(dc$units=="degrees"&dc$rotation=="clock"){
        label.pos=seq(0, (2*4 - 2)*pi/4, length.out = 4)
        labels=(c(0,270,180,90)+conversion.circular(dc$zero,units=dc$units))%%(360)
      }else if(dc$units=="degrees"&dc$rotation=="counter"){
        label.pos=seq(0, (2*4 - 2)*pi/4, length.out = 4)
        labels=(c(0,90,180,270)-conversion.circular(dc$zero,units=dc$units))%%(360)
      }else if(dc$units=="radians"&dc$rotation=="counter"){
        ll=c(0,pi/2,pi,3*pi/2)
        label.pos <- (ll+dc$zero) %% (2*pi)
        labels<-c(0,expression(pi/2),expression(pi),expression(3*pi/2))
      }else if(dc$units=="radians"&dc$rotation=="clock"){
        ll=c(0,pi/2,pi,3*pi/2)
        label.pos <- (-ll + dc$zero) %% (2*pi)
        labels<-c(0,expression(pi/2),expression(pi),expression(3*pi/2))
      }
      
      
    }else{k<-length(labels);label.pos<-seq(0, 2*pi - 2*pi/k, length.out = k)}
    
    
    if(dc$modulo=="pi"){
      
      fx<-dens_x$y
      tx<-dens_x$x
      #if (dc$rotation == "clock") {tx <- -tx}
      tx <- tx + dc$zero
      tx <- tx%%(2 * pi)
      tx <- as.vector(tx)/2
      
      
      r <- sqrt(fx)
      
      or<-order(tx) 
      
      
      radial.plot(
        lengths = c(r[or],r[or]),
        radial.pos=c(tx[or],tx[or]+pi),
        rp.type="p",
        lwd = lwd.dens,
        label.pos <- label.pos,
        radial.lim=c(0,max(r+0.05)),
        show.centroid=FALSE,
        show.grid.labels = FALSE,
        grid.bg="grey95",
        line.col=col.dens,
        labels =labels
      )
      
      
      #if (dc$rotation == "clock") {modes <- -modes}
      modes <- modes + dc$zero
      modes <- modes%%(2 * pi)/2
      testlen <- rep(sqrt(dens_x_modes),2)
      testpos <- c(modes,modes+pi)
      radial.plot(testlen,
                  testpos,
                  add=TRUE,
                  line.col=col.modes,
                  lwd=lwd.modes
      )
      
      
    }else{
      
      fx<-dens_x$y
      tx<-dens_x$x
      
      #if (dc$rotation == "clock") {tx <- -tx}
      #tx <- tx + dc$zero
      tx <- tx%%(2 * pi)
      tx <- as.vector(tx)
      
      
      r <- sqrt(fx)
      
      or<-order(tx) 
      
      radial.plot(
        lengths = r[or],
        radial.pos=tx[or],
        rp.type="p",
        lwd = lwd.dens,
        label.pos=label.pos,
        radial.lim=c(0,max(r+0.05)),
        line.col=col.dens,
        show.centroid=FALSE,
        show.grid.labels = FALSE,
        grid.bg="grey95",
        labels = labels
      )
      
      
      #if (dc$rotation == "clock") {modes <-  -modes}
      #modes <- modes + dc$zero
      modes <- modes%%(2 * pi)
      testlen <- sqrt(dens_x_modes)
      testpos <- modes
      radial.plot(testlen,
                  testpos,
                  add=TRUE,
                  line.col=col.modes,
                  lwd=lwd.modes
      )
      
    }
    
  }
  
  
  return(list(data = data, modes = modes_transformed, bw.modes = bw.modes, bw.dens=kappa))
}

