#' Transition Probability Matrix Utilities
#'
#' @param cdf monotone increasing vector with positive-response probabilities. The number of dose levels $M$ is deduced from vector's length.
#' @param target the design's target response rate.
#' @param repeatNegatives logical, relevant only for k-in-a-row designs: are repeated negative responses needed to level up, or vice versa (repeated positives to level down)? See "Details" for more information.
#' @value for `bcdmat` and `kmatMarg`, an $M\times M$ transition probability matrix. For `pivec`,  an $M$-length vector with stationary/asymptotic visit frequencies.


##################################################################
### Functions that return a matrix

### BCD matrix ##############################

bcdmat<-function(cdf,target)
{
# External validation
validUDinput(cdf,target)

m=length(cdf)
omat=matrix(0,nrow=m,ncol=m)


# Filling in the tridiagonal matrix
# Note: filling off-diagonals already designed to include boundary conditions
# (boundary coded as "diagonal off-diagonal")
if(target<=0.5)
{
	coin=target/(1-target)
# Down moves
	omat[cbind(1:m,c(1,1:(m-1)))]=cdf
# Up moves
	omat[cbind(1:m,c(2:m,m))]=(1-cdf)*coin
# The remainder from 1 goes in the diagonal
	diag(omat)=diag(omat)+(1-cdf)*(1-coin)
	
} else {
	coin=(1-target)/target
# Up moves
	omat[cbind(1:m,c(2:m,m))]=1-cdf
# Down moves
	omat[cbind(1:m,c(1,1:(m-1)))]=cdf*coin
# The remainder from 1 goes in the diagonal
	diag(omat)=diag(omat)+cdf*(1-coin)
}

return(omat)
}

############## K-in-row (geometric) *marginal* stationary matrix 
##############  (one state for each dose)

kmatMarg<-function(cdf,k,repeatNegatives=TRUE)
{
# Finding target from k and direction
kpower=0.5^(1/k)
target=ifelse(repeatNegatives,1-kpower,kpower)
# External validation
validUDinput(cdf,target)

m=length(cdf)
omat=matrix(0,nrow=m,ncol=m)

# useful shorthands
fpower=(1-cdf)^k

# Filling in the tridiagonal matrix
# Note: filling off-diagonals already designed to include boundary conditions
# (boundary coded as "diagonal off-diagonal")
if(target<=0.5)
{
# Down moves
	omat[cbind(1:m,c(1,1:(m-1)))]=cdf
# Up moves
	omat[cbind(1:m,c(2:m,m))]=cdf*fpower/(1-fpower)
# The remainder from 1 goes in the diagonal
	diag(omat)=diag(omat)+1-cdf/(1-fpower)

} else {
# Up moves
	omat[cbind(1:m,c(2:m,m))]=1-cdf
# Down moves
	omat[cbind(1:m,c(1,1:(m-1)))]=(1-cdf)*cdf^k/(1-cdf^k)
# The remainder from 1 goes in the diagonal
	diag(omat)=diag(omat)+1-(1-cdf)/(1-cdf^k)
}

return(omat)
}

############## K-in-row (geometric) *full* stationary matrix 
##############  (with internal states for each dose except top one)

kmatFull<-function(cdf,k,repeatNegatives=TRUE,fluffup=FALSE)
{
# Finding target from k and direction
kpower=0.5^(1/k)
target=ifelse(repeatNegatives,1-kpower,kpower)
# External validation
validUDinput(cdf,target)

m=length(cdf)
mm=(m-1)*k+1
omat=matrix(0,nrow=mm,ncol=mm)

# Note: code already designed to include boundary conditions without explicit exceptions
# Note2: here the diagonal is empty except at boundaries
if(target<=0.5)
{
# Down moves
	omat[cbind(1:mm,rep(c(1,k*(1:(m-1))-k+1),each=k)[1:mm])]=rep(cdf,each=k)[1:mm]
# Up moves
	omat[cbind(1:mm,c(2:mm,mm))]=rep(1-cdf,each=k)[1:mm]

} else {
# Up moves
	omat[cbind(1:mm,c(k+1,rep(k*c(2:(m-1),m-1)+1,each=k)))]=c(1-cdf[1],rep(1-cdf[-1],each=k))
# Down moves
	omat[cbind(1:mm,c(1,1:(mm-1)))]=c(cdf[1],rep(cdf[-1],each=k))
}

return(omat)
}




####################################################
### Functions that return a vector

pivec<-function(cdf,matfun=bcdmat,...)
{
m=length(cdf)
imat=matfun(cdf,...)

vout=cumprod(c(1,imat[cbind(1:(m-1),2:m)]/imat[cbind(2:m,1:(m-1))]))
vout/sum(vout)
}

######

advanceMat<-function(startdose=NULL,cdf,n,matfun=bcdmat,...)
{
require(expm)
m=length(cdf)
## Starting vector
vec0=rep(1/m,m)
if (!is.null(startdose) && startdose %in% 1:m) {
	vec0=rep(0,m)
	vec0[startdose]=1
}
## custom probability vector
if (length(startdose)==m) vec0=startdose/sum(startdose)

vec0 %*% (matfun(cdf,...) %^% (n-1))
}

###

cumulpi<-function(startdose=NULL,cdf,n,matfun=bcdmat,average=TRUE,exclude=1,...)
{
require(expm)
m=length(cdf)
## Starting vector
vec0=rep(1/m,m)
if (!is.null(startdose) && startdose %in% 1:m) {
	vec0=rep(0,m)
	vec0[startdose]=1
}
## custom probability vector
if (length(startdose)==m) vec0=startdose/sum(startdose)

progmat=matfun(cdf,...)

ovec=rep(0,m)
if (is.null(exclude)) exclude=0
for (a in exclude:(n-1)) ovec=ovec+vec0 %*% (progmat %^% a)
if(average) ovec=ovec/(n-exclude)

ovec
}


############################### Auxiliary utilities


validUDinput<-function(cdf,target)
{
if(target<=0 || target>=1) stop("Target has to be in (0,1).\n")
if(min(cdf)<0 || max(cdf)>1 || any(diff(cdf)<0) || var(cdf)==0) stop("cdf should be a CDF.\n")

#ttarg=ifelse(target>0.5,1-target,target)
if(length(cdf)<3) stop ("These designs don't work with <3 dose levels.\n")
}
