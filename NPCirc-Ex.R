pkgname <- "NPCirc"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('NPCirc')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("circsizer.density")
### * circsizer.density

flush(stderr()); flush(stdout())

### Name: circsizer.density
### Title: CircSiZer map for density
### Aliases: circsizer.density
### Keywords: circular density

### ** Examples

# set.seed(2012)
# x <- rcircmix(100,model=7)
# circsizer.density(x, NU=seq(0,50,length=12),type=4,zero=0,clockwise=FALSE)



cleanEx()
nameEx("circsizer.regression")
### * circsizer.regression

flush(stderr()); flush(stdout())

### Name: circsizer.regression
### Title: CircSiZer map for regression
### Aliases: circsizer.regression
### Keywords: circular regression

### ** Examples

# Not run: the code works but it is slow
# set.seed(2012)
# n <- 100
# x <- seq(0,2*pi,length=n)
# y <- sin(x)+sqrt(0.5)*rnorm(n)
# circsizer.regression(x, y, NU=seq(10,60,by=5), title="CircSiZer for regression")



cleanEx()
nameEx("cross.beds1")
### * cross.beds1

flush(stderr()); flush(stdout())

### Name: cross.beds1
### Title: Cross-beds azimuths (I)
### Aliases: cross.beds1
### Keywords: datasets

### ** Examples

data(cross.beds1)



cleanEx()
nameEx("cross.beds2")
### * cross.beds2

flush(stderr()); flush(stdout())

### Name: cross.beds2
### Title: Cross-beds (II)
### Aliases: cross.beds2
### Keywords: datasets

### ** Examples

data(cross.beds2)



cleanEx()
nameEx("cycle.changes")
### * cycle.changes

flush(stderr()); flush(stdout())

### Name: cycle.changes
### Title: Cycle changes
### Aliases: cycle.changes
### Keywords: datasets

### ** Examples

data(cycle.changes)

thaw <- (cycle.changes[,1]==1)
frosting <- (cycle.changes[,1]==-1)

plot(circular(cycle.changes[frosting,2]), shrink=1.08, col=4, stack=TRUE, zero=pi/2,
rotation="clock", axes=FALSE, main="Frosting")
axis.circular(at=circular(seq(0, 7/4*pi, pi/4)), 
labels=c("0h","3h","6h","9h","12h","15h","18h","21h"),
zero=pi/2, rotation="clock")

plot(circular(cycle.changes[thaw,2]), shrink=1.08, col=2, stack=TRUE, zero=pi/2,
rotation="clock", axes=FALSE, main="Thaw")
axis.circular(at=circular(seq(0, 7/4*pi, pi/4)), 
labels=c("0h","3h","6h","9h","12h","15h","18h","21h"),
zero=pi/2, rotation="clock")



cleanEx()
nameEx("dcircmix")
### * dcircmix

flush(stderr()); flush(stdout())

### Name: dcircmix
### Title: Mixtures of circular distributions
### Aliases: dcircmix rcircmix
### Keywords: circular density

### ** Examples

set.seed(2012)
# Linear representation of models M1-M20, each one in a separate window
for (i in 1:20){
dev.new()
f <- function(x) dcircmix(x, model=i)
curve(f, from=0, to=2*pi, main=i)
}

# Circular representation of model vM(0,1)
vM <- function(x) dcircmix(x, model=NULL, dist="vm", param=list(p=1, mu=0, con=1))
curve.circular(vM, n=1000, xlim=c(-1.65,1.65), main="vM(0,1)")
# Random generation from a vM(0,1)
datavM <- rcircmix(50, model=NULL, dist="vm", param=list(p=1, mu=0, con=1))
points(circular(datavM))

# Circular representation of model M18
f18 <- function(x) dcircmix(x, model=18)
curve.circular(f18, n=1000, xlim=c(-1.65,1.65), main="Model 18")
# Random generation from model M8
data18 <- rcircmix(50, model=18)
points(circular(data18))

# Density function and random generation from a mixture of a von Mises and
# a wrapped skew-normal
f <- function(x) dcircmix(x, model=NULL, dist=c("vm","wsn"),
param=list(p=c(0.5,0.5), mu=c(0,pi), con=c(1,1), sk=c(0,10), l=10))
curve.circular(f, n=500, xlim=c(-1.65,1.65))
data <- rcircmix(100, model=NULL, dist=c("vm","wsn"),
param=list(p=c(0.5,0.5), mu=c(0,pi), con=c(1,1), sk=c(0,10)))
points(circular(data))

# Density function and random generation from a mixture of two von Mises and
# two wrapped Cauchy
f <- function(x) dcircmix(x, model=NULL, dist=c("vm","vm","wc","wc"),
param=list(p=c(0.3,0.3,0.2,0.2), mu=c(0,pi,pi/2,3*pi/2), con=c(5,5,0.9,0.9)))
curve.circular(f, n=1000, xlim=c(-1.65,1.65))
data <- rcircmix(100, model=NULL, dist=c("vm","vm","wc","wc"),
param=list(p=c(0.3,0.3,0.2,0.2), mu=c(0,pi,pi/2,3*pi/2), con=c(5,5,0.9,0.9)))
points(circular(data))



cleanEx()
nameEx("dragonfly")
### * dragonfly

flush(stderr()); flush(stdout())

### Name: dragonfly
### Title: Orientations of dragonflies
### Aliases: dragonfly
### Keywords: datasets

### ** Examples

data(dragonfly)
plot(circular(dragonfly), shrink=1.3)
t <- seq(0,2*pi,length=500)
dens <- kern.den.circ(dragonfly$orientation, t, nu=10)
lines(circular(t), dens)



cleanEx()
nameEx("kern.den.circ")
### * kern.den.circ

flush(stderr()); flush(stdout())

### Name: kern.den.circ
### Title: Nonparametric circular kernel density estimation
### Aliases: kern.den.circ
### Keywords: circular density

### ** Examples

##  Estimating the density function of a sample of circular data
set.seed(2012)
n <- 100
x <- rcircmix(n, model=14)
t <- seq(0,2*pi,length=500)
est <- kern.den.circ(x, t, nu=50)
plot(t, dcircmix(t,model=14), ylim=c(0,0.4), type="l", lwd=2, main="Linear representation")
lines(t, est, col=2)
plot(circular(x), shrink=1.3, main="Circular representation")
lines(circular(t), dcircmix(t,model=14), lwd=2)
lines(circular(t), est, col=2)



cleanEx()
nameEx("kern.reg.circ")
### * kern.reg.circ

flush(stderr()); flush(stdout())

### Name: kern.reg.circ
### Title: Nonparametric circular kernel regression estimation
### Aliases: kern.reg.circ
### Keywords: circular regression

### ** Examples

data(speed.wind2)
dir <- rad(speed.wind2$Direction)
vel <- speed.wind2$Speed
nas <- which(is.na(vel))
dir <- dir[-nas]
vel <- vel[-nas]
t <- seq(0,2*pi,length=200)
estLL <- kern.reg.circ(dir, vel, t=t, nu=30)
estNW <- kern.reg.circ(dir, vel, t=t, nu=30, method="NW")

plot(dir, vel, xlab="direction", ylab="speed (m/s)", axes=FALSE)
lines(t, estLL, col=2)
lines(t, estNW, col=3)
axis(1, at=circular(seq(0,2*pi,by=pi/4)),
labels=c("N","NE","E","SE","S","SW","W","NW","N"))
legend("topleft", c("LL","NW"), lty=1, col=2:3)
axis(2)



cleanEx()
nameEx("nu.CV")
### * nu.CV

flush(stderr()); flush(stdout())

### Name: nu.CV
### Title: Cross-validation for density estimation
### Aliases: nu.CV
### Keywords: circular density

### ** Examples

set.seed(2012)
n <- 100
x <- rcircmix(n, model=11)
nu.CV(x, method="LCV", lower=0, upper=20)
nu.CV(x, method="LSCV", lower=0, upper=20)



cleanEx()
nameEx("nu.LSCV.reg")
### * nu.LSCV.reg

flush(stderr()); flush(stdout())

### Name: nu.LSCV.reg
### Title: Least squares cross-validation for circular-linear regression
###   estimation
### Aliases: nu.LSCV.reg
### Keywords: circular regression

### ** Examples

set.seed(2012)
n <- 100
x <- seq(0,2*pi,length=n)
y <- sin(x)+0.2*rnorm(n)
nu.LSCV.reg(x, y, method="LL", lower=1, upper=20)
nu.LSCV.reg(x, y, method="NW", lower=1, upper=20)



cleanEx()
nameEx("nu.boot")
### * nu.boot

flush(stderr()); flush(stdout())

### Name: nu.boot
### Title: Bootstrap method
### Aliases: nu.boot
### Keywords: circular density

### ** Examples

set.seed(2012)
n <- 100
x <- rcircmix(n, model=17)
nu.boot(x, lower=0, upper=20)



cleanEx()
nameEx("nu.pi")
### * nu.pi

flush(stderr()); flush(stdout())

### Name: nu.pi
### Title: Plug-in rule
### Aliases: nu.pi
### Keywords: circular density

### ** Examples

set.seed(2012)
n <- 100
x <- rcircmix(n,model=18)
nu.pi(x, M=3)
nu.pi(x, outM=TRUE)  # Using AIC



cleanEx()
nameEx("nu.rt")
### * nu.rt

flush(stderr()); flush(stdout())

### Name: nu.rt
### Title: Rule of thumb
### Aliases: nu.rt
### Keywords: circular density

### ** Examples

set.seed(2012)
n <- 100
x <- rcircmix(n,model=7)
nu.rt(x) 
nu.rt(x, robust=TRUE)



cleanEx()
nameEx("speed.wind")
### * speed.wind

flush(stderr()); flush(stdout())

### Name: speed.wind
### Title: Wind speed and wind direction data
### Aliases: speed.wind speed.wind2
### Keywords: datasets

### ** Examples

data(speed.wind2)

t <- seq(0,2*pi,length=500)
dir <- rad(speed.wind2$Direction)

# Density
plot(circular(dir), stack=TRUE, zero=pi/2, rotation="clock", axes=FALSE)
axis.circular(at=circular(seq(0,2*pi,by=pi/4)),
labels=c("N","NE","E","SE","S","SW","W","NW","N"), zero=pi/2, rotation="clock")
rose.diag(circular(dir), bins=16, add=TRUE, axes=FALSE, zero=pi/2, rotation="clock")

rose.diag(circular(dir), bins=16, shrink=1.1, axes=FALSE, zero=pi/2, rotation="clock")
axis.circular(at=circular(seq(0,7*pi/4,by=pi/4)),
labels=c("N","NE","E","SE","S","SW","W","NW"), zero=pi/2, rotation="clock")
lines(circular(t), kern.den.circ(dir,t,nu=1), lwd=2, lty=2, zero=pi/2, rotation="clock")
lines(circular(t), kern.den.circ(dir,t,nu=10), lwd=2, lty=1, zero=pi/2, rotation="clock")
lines(circular(t), kern.den.circ(dir,t,nu=60), lwd=2, lty=3, zero=pi/2, rotation="clock")

# Regression
vel <- speed.wind2$Speed
nas <- which(is.na(vel))
dir <- dir[-nas]
vel <- vel[-nas]
radial.plot(vel, dir, rp.type="s", start=pi/2, clockwise=TRUE, point.col="gray",
radial.lim=c(0,15), label.pos=seq(0,7*pi/4,by=pi/4),
labels=c("N","NE","E","SE","S","SW","W","NW"))
radial.plot(as.numeric(kern.reg.circ(dir,vel,t,nu=1,method="LL")), t, rp.type="p",
add=TRUE, start=pi/2, clockwise=TRUE, radial.lim=c(0,15), lwd=2, lty=2)
radial.plot(as.numeric(kern.reg.circ(dir,vel,t,nu=10,method="LL")),t, rp.type="p",
add=TRUE, start=pi/2, clockwise=TRUE, radial.lim=c(0,15), lwd=2, lty=1)
radial.plot(as.numeric(kern.reg.circ(dir,vel,t,nu=60,method="LL")),t, rp.type="p",
add=TRUE, start=pi/2, clockwise=TRUE, radial.lim=c(0,15), lwd=2, lty=3)




cleanEx()
nameEx("temp.wind")
### * temp.wind

flush(stderr()); flush(stdout())

### Name: temp.wind
### Title: Temperature and wind direction data
### Aliases: temp.wind
### Keywords: datasets

### ** Examples

data(temp.wind)

winddir <- rad(temp.wind[,4]) # wind direction in radians
temp <- temp.wind[,3]
nas <- which(is.na(winddir))
winddir <- winddir[-nas] 
temp <- temp[-nas]

# value of the smoothing parameter selected using the function nu.LSCV.reg 
# with method="LL"
nu <- 3.41 
t <- seq(0,2*pi,length=100)
est <- kern.reg.circ(winddir, temp, t, nu=nu, method="LL")

# Circular representation
radial.plot(temp, winddir, rp.type="s", labels=c("N","NE","E","SE","S","SW","W","NW"),
start=pi/2, clockwise=TRUE, label.pos=seq(0,7*pi/4,by=pi/4), lwd=2, point.col="grey50",
radial.lim=c(-10,15))
radial.plot(as.vector(est), rp.type="p", start=pi/2, clockwise=TRUE, lwd=2,
radial.lim=c(-10,15), add=TRUE)

# Linear representation
plot(t, est, type="l", xlab="", ylab="Temperatute (ºC)", axes=FALSE)
axis(1, labels=c("N","NE","E","SE","S","SW","W","NW","N"), at=seq(0,2*pi,by=pi/4))
axis(2)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
