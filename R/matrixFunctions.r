#' Transition Probability Matrices for Up-and-Down Designs
#' 
#' 
#' Transition Probability Matrices for Common Up-and-Down Designs
#'
#' @details
#' Up-and-Down designs (UDDs) generate random walk behavior, whose theoretical properties can be summarized via a transition probability matrix (TPM). Given the number of doses \eqn{M}, and the value of the cdf \eqn{F} (i.e., the positive-response probabilities) at each dose, the specific UDD rules uniquely determine the TPM.
#' 
#' The utilities described here calculate the TPMs for the most common and simplest UDDs:
#' 
#'  - The k-in-a-row or ``fixed staircase`` design common in sensory studies: `kmatMarg(), kmatFull()` (see Note). The design parameters are k, a natural number, and whether k negative responses are required for dose transition, or k positive responses. The former is for targets below the median and vice versa.
#'  - The Durham-Flournoy Biased Coin Design: `bcdmat()` (. This design can target any percentile via the `target` argument.
#'  - The original *"classical"* median-targeting UDD: `cudmat()` (Dixon and Mood, 1948). This is simply a wrapper for `bcdmat()` with `target` set to 0.5.
#'  - Cohort or group UDD: `gudmat()` (Gezmu and Flournoy, 2006).
#'  
#'  
#'  @notes As Gezmu (1996) discovered and Oron and Hoff (2009) further extended, k-in-a-row UDDs with \eqn{k>1} generate a random walk with internal states. Their full TPM is therefore larger than \eqn{M\times M.} However, in terms of random-walk behavior, most salient properties are better represented via an \eqn{M\times M} matrix analogous to those of the other designs, with transition probabilities marginalized over internal states using their asymptotic frequencies. This matrix is provided by `kmatMarg()`, while `kmatFull()` returns the full matrix including internal states.
#'  
#'  Also, in `kmatFull()` there are two matrix-size options. Near one of the boundaries (upper boundary with `repeatNegatives = TRUE` and vice versa), the most extreme \eqn{k} internal states are practically indistinguishable so in some sense they don't really exist. Using the `fluffup` argument, user can choose between having a more pleasantly symmetric (but a bit misleading) full \eqn{Mk\times Mk} matrix, or reducing it to its practically true size by \eqn{k-1} rows and columns.
#'  
#'
#' @param cdf monotone increasing vector with positive-response probabilities. The number of dose levels $M$ is deduced from vector's length.
#' @param target the design's target response rate (`bcdmat()` only).
#' @param k the number of consecutive identical responses required for dose transitions (k-in-a-row functions only).
#' @param repeatNegatives logical: are k repeated negative responses needed to level up, or vice versa (k repeated positives to level down)? k-in-a-row functions only. See "Details" for more information.
#' @param fluffup logical (`kmatFull` only): in the full k-in-a-row internal-state representation, should we *"fluff"* the matrix up so that it has \eqn{Mk} rows and columns (`TRUE`, default), or exclude \eqn{k-1} "phantom" states near one of the boundaries?
#' @param cohort,lower,upper (`gudmat` only): the cohort (group) size, how many positive responses are allowed for a move upward, and how many are required for a move downward, respectively. For example `cohort=3, lower=0, upper=2` evaluates groups of 3 observations at a time, moves up if none are positive, down if \eqn{>=2} are positive, and repeats the same dose with 1 positive.

#' @return An \eqn{M\times M} transition probability matrix, except for `kmatFull()` with \eqn{k>1} which returns a larger square matrix. 

#' @seealso 
#'  - \code{\link{k2targ}}, \code{\link{targ2k}} to find the k-in-a-r-w target-response rate for specific k and vice versa
#'  - \code{\link{gudtarg}} to calculate group UDD target-response rate given a `cohort,lower,upper` combination


#' @references 
#'  - Dixon WJ, Mood AM. A method for obtaining and analyzing sensitivity data. *J Am Stat Assoc.* 1948;43:109-126.
#'  - Durham SD, Flournoy N. Random walks for quantile estimation. In: *Statistical Decision Theory and Related Topics V* (West Lafayette, IN, 1992). Springer; 1994:467-476.
#'  - Gezmu M. The Geometric Up-and-Down Design for Allocating Dosage Levels. PhD Thesis. American University; 1996.
#'  - Gezmu M, Flournoy N. Group up-and-down designs for dose-finding. *J Stat Plan Inference.* 2006;136(6):1749-1764.
#'  - Oron AP, Hoff PD. The k-in-a-row up-and-down design, revisited. *Stat Med.* 2009;28:1805-1820.
#'  
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' @export

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
downmove=cdf
omat[cbind(1:m,c(1,1:(m-1)))]=downmove
# Up moves
upmove=(1-cdf)*coin
omat[cbind(1:m,c(2:m,m))]=upmove
# The remainder from 1 goes in the diagonal
diag(omat)=diag(omat)+1-downmove-upmove
	
} else {
	coin=(1-target)/target
# Up moves
upmove=1-cdf
omat[cbind(1:m,c(2:m,m))]=upmove
# Down moves
downmove=cdf*coin
omat[cbind(1:m,c(1,1:(m-1)))]=downmove
# The remainder from 1 goes in the diagonal
diag(omat)=diag(omat)+1-downmove-upmove
}

return(omat)
}

# Classical, using bcdmat()

#' @rdname bcdmat
#' @export

cudmat <- function(cdf) bcdmat(cdf, target = 1/2)


############## K-in-row (geometric) *marginal* stationary matrix 
##############  (one state for each dose)

#' @rdname bcdmat
#' @export

kmatMarg <- function(cdf,k,repeatNegatives=TRUE)
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
downmove=cdf
omat[cbind(1:m,c(1,1:(m-1)))]=downmove
# Up moves
upmove=cdf*fpower/(1-fpower)
omat[cbind(1:m,c(2:m,m))]=upmove
# The remainder from 1 goes in the diagonal
diag(omat)=diag(omat)+1-downmove-upmove

} else {

# Up moves
upmove=1-cdf
omat[cbind(1:m,c(2:m,m))]=upmove
# Down moves
downmove=(1-cdf)*cdf^k/(1-cdf^k)
omat[cbind(1:m,c(1,1:(m-1)))]=downmove
# The remainder from 1 goes in the diagonal
diag(omat)=diag(omat)+1-downmove-upmove
}

return(omat)
}


############## K-in-row (geometric) *full* stationary matrix 
##############  (with internal states for each dose except top one)

#' @rdname bcdmat
#' @export

kmatFull<-function(cdf,k,repeatNegatives=TRUE,fluffup=FALSE)
{
# Finding target from k and direction
kpower=0.5^(1/k)
target=ifelse(repeatNegatives,1-kpower,kpower)
# External validation
validUDinput(cdf, target)

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

### GUD matrix ##############################

#' @rdname bcdmat
#' @export

gudmat<-function(cdf,cohort,lower,upper)
{
## Validation (lots!)
if(min(cdf)<0 || max(cdf)>1 || any(diff(cdf)<0) || var(cdf)==0) stop("cdf should be a CDF.\n")
if(length(cdf)<3) stop ("These designs don't work with <3 dose levels.\n")
if(cohort!=round(cohort) || upper!=round(upper) || lower!=round(lower)) stop('Design parameters must be integers.\n')
if(cohort<1 || upper<1 || cohort<upper) stop('upper<=cohort and both must be positive.\n')
if(lower<0 || upper<=lower) stop('lower<upper and lower cannot be negative.\n')
# /validation

m=length(cdf)
omat=matrix(0,nrow=m,ncol=m)

## Filling in the tridiagonal matrix
# Note: filling off-diagonals already designed to include boundary conditions
# (boundary coded as "diagonal off-diagonal")

# Down moves
downmove=pbinom(upper-1,size=cohort,prob=cdf,lower.tail=FALSE)
omat[cbind(1:m,c(1,1:(m-1)))]=downmove
# Up moves
upmove=pbinom(lower,size=cohort,prob=cdf)
omat[cbind(1:m,c(2:m,m))]=upmove
# The remainder from 1 goes in the diagonal
diag(omat)=diag(omat)+1-downmove-upmove
	
return(omat)
}

####################################################
### Functions that return a vector
### They go in a different help file

pivec<-function(cdf,matfun=bcdmat,...)
{
m=length(cdf)
imat=matfun(cdf,...)

vout=cumprod(c(1,imat[cbind(1:(m-1),2:m)]/imat[cbind(2:m,1:(m-1))]))
vout/sum(vout)
}

######

advanceVec<-function(startdose=NULL,cdf,n,designMat=bcdmat,...)
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

vec0 %*% (designMat(cdf,...) %^% (n-1))
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




# For `pivec`,  an $M$-length vector with stationary/asymptotic visit frequencies.

