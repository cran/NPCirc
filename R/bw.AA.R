
bw.AA=function(x,deriv.order=0,method = c("ste","dpi"),nstage=2,kernel="vonmises",M=NULL,commonkappa=TRUE,Q1=NULL,Q2=NULL,lower=NULL,upper=NULL,tol=NULL,approximate=NULL){

  data=conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL

  if (!is.numeric(data)) stop("Argument 'data' must be numeric")
  if (sum(is.na(data)) > 0) warning("Missing values were removed")
  data=data[!is.na(data)]

  n=length(data)

  if (n<2) stop("At least 2 data points are needed")


  if (!is.null(deriv.order)) {
    if (deriv.order!=as.integer(deriv.order)|deriv.order<0|length(deriv.order)!=1) {
      warning("The provided value of 'deriv.order' is not valid, 'deriv.order' must be a positive integer. Default value of 'deriv.order' was used")
      deriv.order=0
    }
  }

  method=match.arg(method)

  if (!is.null(nstage)) {
    if (nstage!=as.integer(nstage)|nstage<0|length(nstage)!=1) {
      warning("The provided value of 'nstage' is not valid, 'nstage' must be a non-negative integer. Default value of 'nstage' was used")
      nstage=2
    }
  }


  # Check if the employed kernel is available at cosmoments
  kernel.check=cosmoments(0.1,1,kernel)

  if (!is.null(M)) {
    if (sum(M!=as.integer(M))>0|sum(M<1)>0) {
      warning("The provided value of 'M' is not valid, 'M' must be a positive integer. 'M' was selected by AIC")
      M = NULL
    }
  }

  if(is.null(M)){
    M=1:5
  }

  if (is.null(M) & method == "ste") {
    M=1
  }

  if (is.null(M) & method == "dpi") {
    M=1:5
  }

  if (commonkappa != T & commonkappa != F) {
    warning("Argument 'commonkappa' must be T or F. Default value of 'commonkappa' was used")
    commonkappa=T
  }



  if(!is.null(Q1)){
    if (!is.numeric(Q1)) {
      warning("argument 'Q1' must be numeric. Default value was used")
      Q1=NULL
    }}

  if(!is.null(Q2)){
    if (!is.numeric(Q2)) {
      warning("argument 'Q2' must be numeric. Default value was used")
      Q1=NULL
    }}

  seqr=2*(deriv.order+1:nstage)+2

  if(!is.null(Q1)){
    if(sum(is.na(Q1[seqr]))>0){
      warning("There are missing needed elements on the vector Q1. Default values were used")
      Q1=NULL
    }
  }

  if(is.null(Q1)){
    Q1=numeric()
    # Q1 (the only needed elements are the ones in the positions seqr)
    if(kernel=="vonmises"|kernel=="wrappednormal"){
      Q1[seqr]=(-1)^(seqr/2)*factorial(seqr)/(2^(seqr/2)*factorial(seqr/2)*sqrt(2*pi))
    }else{
      stop("vector 'Q1' needs to be imposed for the employed kernel")
    }

  }

  if(is.null(Q2)){
    sr=deriv.order
    if(kernel=="vonmises"|kernel=="wrappednormal"){
      Q2=factorial(2*sr)/(2^(2*sr+1)*factorial(sr)*sqrt(pi))
    }else{
      stop("argument 'Q2' needs to be imposed for the employed kernel")
    }
  }

  if (!is.null(approximate) & method == "ste") {
    warning("Argument 'approximate' is not employed for this method")
  }


  if (is.null(approximate) & method == "dpi") {
    approximate=T
  }


  if (!is.null(lower) & method == "dpi") {
    warning("Argument 'lower' is not employed for this method")
  }

  if (is.null(lower) & method == "ste") {
    lower=10^(-10)
    if(kernel=="vonmises"){lower=10^(-3)}
  }

  if (!is.null(upper) & method == "dpi") {
    warning("Argument 'upper' is not employed for this method")
  }

  if (is.null(upper) & method == "ste") {
    upper=0.5
    if(kernel=="vonmises"){upper=pi^2/3}
  }

  if (!is.null(tol) & method == "dpi") {
    warning("Argument 'tol' is not employed for this method")
  }

  if (is.null(tol) & method == "ste") {
    tol=.Machine$double.eps^0.25
  }


  if(method == "dpi"){
    return(bw.dpi(data,deriv.order,nstage,kernel,M,commonkappa,approximate,Q1,Q2))
  }

  if(method == "ste"){
    return(bw.ste(data,deriv.order,nstage,kernel,M,commonkappa,lower,upper,tol,Q1,Q2))
  }

}
