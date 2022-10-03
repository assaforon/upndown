#' Transition Probability Matrices for Up-and-Down Designs
#' 
#' Transition Probability Matrices for Common Up-and-Down Designs
#'
#' @details
#' Up-and-Down designs (UDDs) generate random walk behavior, whose theoretical properties can be summarized via a transition probability matrix (TPM). Given the number of doses \eqn{M}, and the value of the cdf \eqn{F} at each dose (i.e., the positive-response probabilities), the specific UDD rules uniquely determine the TPM.
#' 
#' The utilities described here calculate the TPMs of the most common and simplest UDDs:
#' 
#'  - The k-in-a-row or **fixed staircase** design common in sensory studies: `kmatMarg(), kmatFull()` (Gezmu, 1996; Oron and Hoff, 2009; see Note). Design parameters are k, a natural number, and whether k negative responses are required for dose transition, or k positive responses. The former is for targets below the median and vice versa.
#'  - The Durham-Flournoy Biased Coin Design: `bcdmat()`. This design can target any percentile via the `target` argument (Durham and Flournoy, 1994).
#'  - The original *"classical"* median-targeting UDD: `classicmat()` (Dixon and Mood, 1948). This is simply a wrapper for `bcdmat()` with `target` set to 0.5.
#'  - Cohort or group UDD: `gudmat()`, with three design parameters for the group size and the up/down rule thresholds  (Gezmu and Flournoy, 2006).
#'  
#'  
#' @note 
#' As Gezmu (1996) discovered and Oron and Hoff (2009) further extended, k-in-a-row UDDs with \eqn{k>1} generate a random walk *with internal states*. Their full TPM is therefore larger than \eqn{M\times M.} However, in terms of random-walk behavior, most salient properties are better represented via an \eqn{M\times M} matrix analogous to those of the other designs, with transition probabilities marginalized over internal states using their asymptotic frequencies. This matrix is provided by `kmatMarg()`, while `kmatFull()` returns the full matrix including internal states.
#'  
#'  Also, in `kmatFull()` there are two matrix-size options. Near one of the boundaries (upper boundary with `lowTarget = TRUE`, and vice versa), the most extreme \eqn{k} internal states are practically indistinguishable, so in some sense only one of them really exists. Using the `fluffup` argument, users can choose between having a more aesthetically symmetric (but a bit misleading) full \eqn{Mk\times Mk} matrix, or reducing it to its effectivelly true size by removing \eqn{k-1} rows and columns.
#'  
#'
#' @param cdf monotone increasing vector with positive-response probabilities. The number of dose levels $M$ is deduced from vector's length.
#' @param target the design's target response rate (`bcdmat()` only).
#' @param k the number of consecutive identical responses required for dose transitions (k-in-a-row functions only).
#' @param lowTarget logical k-in-a-row functions only: is the design targeting below-median percentiles, with \eqn{k} repeated negative responses needed to level up and only one to level down - or vice versa? Default `FALSE`. See "Details" for more information.
#' @param fluffup logical (`kmatFull` only): in the full k-in-a-row internal-state representation, should we *"fluff"* the matrix up so that it has \eqn{Mk} rows and columns (`TRUE`, default), or exclude \eqn{k-1} "phantom" states near one of the boundaries?
#' @param cohort,lower,upper `gudmat` only: the cohort (group) size, how many positive responses are allowed for a move upward, and how many are required for a move downward, respectively. For example `cohort=3, lower=0, upper=2` evaluates groups of 3 observations at a time, moves up if none are positive, down if \eqn{>=2} are positive, and repeats the same dose with 1 positive.

#' @return An \eqn{M\times M} transition probability matrix, except for `kmatFull()` with \eqn{k>1} which returns a larger square matrix. 

#' @seealso 
#'  - \code{\link{k2targ}}, \code{\link{targ2k}} to find the k-in-a-row target-response rate for specific k and vice versa
#'  - \code{\link{gudtarg}} to calculate group UDD target-response rate given a `cohort,lower,upper` combination
#'  - \code{\link{pivec}}, \code{\link{advancevec}}, \code{\link{cumulpi}}, which provide probability vectors of dose-allocation distributions using Up-and-Down TPMs.


#' @references 
#'  - Dixon WJ, Mood AM. A method for obtaining and analyzing sensitivity data. *J Am Stat Assoc.* 1948;43:109-126.
#'  - Durham SD, Flournoy N. Random walks for quantile estimation. In: *Statistical Decision Theory and Related Topics V* (West Lafayette, IN, 1992). Springer; 1994:467-476.
#'  - Gezmu M. The Geometric Up-and-Down Design for Allocating Dosage Levels. PhD Thesis. American University; 1996.
#'  - Gezmu M, Flournoy N. Group up-and-down designs for dose-finding. *J Stat Plan Inference.* 2006;136(6):1749-1764.
#'  - Oron AP, Hoff PD. The k-in-a-row up-and-down design, revisited. *Stat Med.* 2009;28:1805-1820.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50.
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

##### Classical, using bcdmat()

#' @rdname bcdmat
#' @export

classicmat <- function(cdf) bcdmat(cdf, target = 1/2)


############## K-in-row (geometric) *marginal* stationary matrix 
##############  (one state for each dose)

#' @rdname bcdmat
#' @export

kmatMarg <- function(cdf, k, lowTarget)
{
#### Validation and prep

checkCDF(cdf)
checkNatural(k, parname = 'k', toolarge = 30)  
# Finding target from k and direction
kpower=0.5^(1/k)
target=ifelse(lowTarget,1-kpower,kpower)
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
##############  (with internal states)

#' @rdname bcdmat
#' @export

kmatFull<-function(cdf, k, lowTarget, fluffup = FALSE)
{
#### Validation and prep

checkCDF(cdf)
checkNatural(k, parname = 'k', toolarge = 30)  
# Finding target from k and direction
kpower=0.5^(1/k)
target=ifelse(lowTarget,1-kpower,kpower)
# External validation
validUDinput(cdf, target)
# A bit more validation/prep:
if(k==1) fluffup = FALSE # fluffup is irrelevant w/k=1

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
	
# Optionally dding "semi-dummy" rows+columns 
#   to get a "nicer" m*k square matrix

	if(fluffup)
	{
	  mm = m * k
	  d = nrow(omat) - 1 # shorthand for the m-1 levels we won't touch
	  omat = rbind(omat[1:d, 1:d], matrix(0, nrow = k, ncol = d))
	  omat = cbind(omat, matrix(0, nrow = d+k, ncol = k))
# Down moves
	  omat[ cbind( (d+1):mm, d-k+1 ) ] = cdf[m]
# "Up" moves (really, meaningless internal-state increments)
# The first one got deleted in the expansion
	  omat[d, d+1] = 1 - cdf[m-1]
	  omat[ cbind( (d+1):mm,c((d+2):mm,mm) ) ] = 1 - cdf[m]
	}                           # end fluffup; no 'else' for this one

} else {       #  target > 0.5

# Up moves
	omat[cbind(1:mm,c(k+1,rep(k*c(2:(m-1),m-1)+1,each=k)))]=c(1-cdf[1],rep(1-cdf[-1],each=k))
# Down moves
	omat[cbind(1:mm,c(1,1:(mm-1)))]=c(cdf[1],rep(cdf[-1],each=k))

	if(fluffup)
	{
	  mm = m * k
	  d = nrow(omat) - 1 # shorthand for the m-1 levels we won't touch
	  omat = rbind( matrix(0, nrow = k, ncol = d), omat[2:(d+1), 2:(d+1)] )
	  omat = cbind(matrix(0, nrow = d+k, ncol = k), omat)
	  # Up moves
	  omat[ cbind( 1:k, k+1 ) ] = 1 -cdf[1]
	  # "Down" moves (really, meaningless internal-state decrements)
	  # The first one got deleted in the expansion
	  omat[k+1, k] = cdf[2]
	  omat[ cbind( 1:k, c( 1,1:(k-1) ) )  ] = cdf[1]
	}              
}

return(omat)
}

### GUD matrix ##############################

#' @rdname bcdmat
#' @export

gudmat<-function(cdf, cohort, lower, upper)
{
## Validation (lots!)
checkCDF(cdf)
checkNatural(c(cohort, lower+1, upper), parname = 'cohort, lower+1, upper', toolarge = 50)  
if(cohort<upper || upper<=lower) stop('Order must be lower < upper <= cohort.\n')
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

