

kern.den.circ=function(x,z=NULL,bw="AA",deriv.order=0,kernel="vonmises",na.rm = FALSE, from = circular(0), to = circular(2 * pi),n = 512,control.circular=list()){

  # Inherit from density.circular function (for compatibility)
  data<-x
  name<-deparse(substitute(x))
  if (!is.numeric(from))
    stop("argument 'from' must be numeric")
  if (!is.numeric(to))
    stop("argument 'to' must be numeric")
  if (!is.finite(from))
    stop("non-finite `from'")
  if (!is.finite(to))
    stop("non-finite `to'")
  if (!is.numeric(n))
    stop("argument 'n' must be numeric")
  n <- round(n)
  if (n <= 0)
    stop("argument 'n' must be integer and positive")
  if (!is.numeric(x))
    stop("argument 'x' must be numeric")
  if (!is.null(z) && is.circular(z)) {
    datacircularp <- circularp(z)
  }
  else if (is.circular(x))
    datacircularp <- circularp(x)
  else {
    datacircularp <- list(type = "angles", units = "radians",
                          template = "none", modulo = "asis", zero = 0,
                          rotation = "counter")
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
  if (dc$modulo == "pi")
    stop("The function does not work yet for modulo='pi'")
  x <- conversion.circular(x, units = "radians", zero = 0,
                           rotation = "counter")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  from <- conversion.circular(from, units = "radians",
                              zero = 0, rotation = "counter")
  attr(from, "class") <- attr(from, "circularp") <- NULL
  to <- conversion.circular(to, units = "radians", zero = 0,
                            rotation = "counter")
  attr(to, "class") <- attr(to, "circularp") <- NULL
  x <- as.vector(x)
  x.na <- is.na(x)
  if (any(x.na)) {
    if (na.rm)
      x <- x[!x.na]
    else stop("x contains missing values")
  }
  x.finite <- is.finite(x)
  if (any(!x.finite)) {
    x <- x[x.finite]
  }
  nx <- length(x)
  if (is.null(z)) {
    z <- circular(seq(from = from, to = to, length = n))
  }
  else {
    if (!is.numeric(z))
      stop("argument 'z' must be numeric")
    namez <- deparse(substitute(z))
    z.na <- is.na(z)
    if (any(z.na)) {
      if (na.rm) {
        z <- z[!z.na]
      }
      else {
        stop("z contains missing values")
      }
    }
    z.finite <- is.finite(z)
    if (any(!z.finite)) {
      z <- z[z.finite]
    }
  }
  zz <- conversion.circular(z, dc$units, dc$type, dc$template,
                            dc$modulo, dc$zero, dc$rotation)
  z <- conversion.circular(z, units = "radians", zero = 0,
                           rotation = "counter")
  attr(z, "class") <- attr(z, "circularp") <- NULL
  z <- as.vector(z)



  if (!is.null(deriv.order)) {
    if (deriv.order!=as.integer(deriv.order)|deriv.order<0|length(deriv.order)!=1) {
      warning("The provided value of 'deriv.order' is not valid, 'deriv.order' must be a positive integer. Default value of 'deriv.order' was used")
      deriv.order=0
    }
  }


  # Check if the employed kernel is available at cosmoments
  kernel.check=cosmoments(0.1,1,kernel)




  if(is.character(bw)){
    if((kernel=="vonmises")&(deriv.order==0)){
      bw <- switch(tolower(bw), aa = bw.AA(data), ste = bw.AA(data,method="ste"), dpi = bw.AA(data,method="dpi"), rt = bw.rt(data), pi = bw.pi(data), cv = bw.CV(data), boot = bw.boot(data), lcv = bw.CV(data), lscv = bw.CV(data,method="LSCV"),  warning("Unknown bandwidth rule. Default value of 'bw' was used"))
    }else{

      bw <- switch(tolower(bw), aa = bw.AA(data,deriv.order=deriv.order,kernel=kernel), ste = bw.AA(data,deriv.order=deriv.order,method="ste",kernel=kernel), dpi = bw.AA(data,deriv.order=deriv.order,method="dpi",kernel=kernel), rt = bw.AA(data,deriv.order=deriv.order,method="dpi",nstage=0,kernel=kernel,M=1), pi = bw.AA(data,deriv.order=deriv.order,method="dpi",nstage=0,kernel=kernel), warning("The bandwidth rule is unknown. Default value of 'bw' was used"))

    }
  }

  if(is.character(bw)){
    bw=bw.AA(data,deriv.order=deriv.order,kernel=kernel)
  }

  if (!is.numeric(bw))
    stop("argument 'bw' and 'adjust' must be numeric")
  if (!is.finite(bw))
    stop("non-finite `bw'")
  if (bw < 0)
    stop("`bw' is not positive.")



  if(kernel=="vonmises"){
    matkdde=ddvonmises(z,data,bw,deriv.order)
  }else{
    matkdde=dcirculardr(z,data,bw,deriv.order,kernel)
  }

  y=rowSums(matkdde)/length(data)

  structure(list(data = data, x = zz, y = y, bw = bw, n = nx,
                 kernel = kernel, call = match.call(), data.name = name,
                 has.na = FALSE), class = "density.circular")

}




