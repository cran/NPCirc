check_pch <- function(pch, arg_name = "pch") {
  # Must be length 1
  if (length(pch) != 1) {
    stop(sprintf("'%s' must be a single value", arg_name))
  }
  
  # Numeric case
  if (is.numeric(pch)) {
    if (pch < 0 || pch > 25) {
      stop(sprintf("'%s' numeric value must be between 0 and 25", arg_name))
    }
    return(as.numeric(pch))
  }
  
  # Character case
  if (is.character(pch)) {
    if (nchar(pch) != 1) {
      stop(sprintf("'%s' character must be a single character", arg_name))
    }
    return(pch)
  }
  
  # Other types are invalid
  stop(sprintf("'%s' must be numeric (0-25) or a single character", arg_name))
}

check_lwd <- function(lwd, arg_name = "lwd") {
  # Must be numeric
  if (!is.numeric(lwd)) {
    stop(sprintf("'%s' must be numeric", arg_name))
  }
  
  # Must be length 1
  if (length(lwd) != 1) {
    stop(sprintf("'%s' must be a single value", arg_name))
  }
  
  # Must be positive
  if (lwd <= 0) {
    stop(sprintf("'%s' must be positive", arg_name))
  }
  
  # Optionally, convert integer to numeric
  lwd <- as.numeric(lwd)
  
  return(lwd)
}

segments.circular<-function (x, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, shrink = 1, 
                             plot.info = NULL, zero = NULL, rotation = NULL, ...) 
{
  if (na.rm) 
    x <- x[!is.na(x)]
  if (length(x) == 0) {
    warning("No observations (at least after removing missing values)")
    return(NULL)
  }
  xcircularp <- attr(as.circular(x), "circularp")
  if (is.null(plot.info)) {
    if (is.null(zero)) 
      zero <- xcircularp$zero
    if (is.null(rotation)) 
      rotation <- xcircularp$rotation
  }
  else {
    zero <- plot.info$zero
    rotation <- plot.info$rotation
  }
  x <- conversion.circular(x, units = "radians")
  attr(x, "class") <- attr(x, "circularp") <- NULL
  if (rotation == "clock") 
    x <- -x
  x <- x + zero
  x <- x%%(2 * pi)
  x <- as.vector(x)
  if (is.null(y)) 
    y <- rep(1, length(x))
  y <- as.vector(y)
  if (length(y) != length(x)) 
    stop("'y' must have the same length of 'x'")
  y <- y * shrink
  if (length(x0) != length(x)) 
    x0 <- rep(x0, length(x))
  if (length(y0) != length(x)) 
    y0 <- rep(y0, length(x))
  x1 <- x0 + y * cos(x)
  y1 <- y0 + y * sin(x)
  segments(x0, y0, x1, y1, ...)
}


circ.meanshift<-function(phi,kappa,phi0,tol){
  
  res<-1000
  cont<-0
  
  while(abs(res)>tol){
    phil<-atan2(sum(exp(kappa*cos(phi0-phi))*sin(phi)),sum(exp(kappa*cos(phi0-phi))*cos(phi)))
    res<-1-cos(phil-phi0)
    phi0<-phil
    cont<-cont+1
  }
  
  return(phil)
}


##Actually performing mode estimation:

mode_estimation<-function(x,bw=NULL,tol=0.00000000001){
  
  if(is.null(bw)){kappa<-bw.AA(x,deriv.order=1)}else{kappa<-bw}
  modes<-numeric(length(x))
  for (i in 1:length(x)){
    modes[i]<-circ.meanshift(x,kappa,x[i],tol)
  }
  z<-round(modes,2)
  modes_unique<-unique(round(modes,2)) # final estimated modes
  
  k<-length(modes_unique)
  colors<-numeric(k)
  
  base_colors <- c(3, 6, 7, 8, 9, 2, 4, "forestgreen", 5)
  
  # If there are more modes than the predefined colors, add extra colors
  if(k > length(base_colors)) {
    # Use a palette of R colors, excluding the ones already in base_colors
    extra_colors <- setdiff(rainbow(k), base_colors)  # or other palette
    all_colors <- c(base_colors, extra_colors[1:(k - length(base_colors))])
  } else {
    all_colors <- base_colors[1:k]
  }
  
  # Assign colors according to modes_unique
  for(i in 1:k) {
    colors[z == modes_unique[i]] <- all_colors[i]
  }
  
  colors_numeric <- as.numeric(factor(colors, levels = unique(colors)))
  
  return(list(modes=modes_unique,bw=kappa,clust=colors_numeric,cols=colors))
  
}



