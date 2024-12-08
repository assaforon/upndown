#' Key Probability Vectors of Up-and-Down Designs
#'
#' Dose-allocation probability vectors that quantify the instantaneous, cumulative, and asymptotic behavior of Up-and-Down designs.
#'
#' @details
#' Up-and-Down designs (UDDs) generate random walk behavior, which concentrates doses around the target quantile. Asymptotically, dose allocations follow a stationary distribution \eqn{\boldsymbol{\pi}} which can be calculated given the number of doses \eqn{M}, and the value of the cdf \eqn{F} at each dose (i.e., the positive-response probabilities), and the specific UDD rules. No matter the starting dose, the allocation distribution converges to \eqn{\boldsymbol{\pi}} at a geometric rate (Diaconis and Stroock, 1991).  
#' 
#' Three functions are offered:
#' 
#'  - `pivec()` returns \eqn{\boldsymbol{\pi}}.
#'  - `currentvec()` returns the current (instantaneous) allocation distribution at step `n`, using the formula from Diaconis and Stroock (1991). 
#'  - `cumulvec()` returns the *cumulative* allocations, i.e., the expected proportions (or counts) of allocations during the experiment after `n` observations. This function is perhaps of greatest practical use. 
#'  
#' All functions first calculate the transition probability matrix (TPM), by calling one of the functions described under \code{\link{bcdmat}}. See that help page for more details.
#' 
#' @return A vector of allocation frequencies/probabilities for the doses, summing up to 1. Exception: `cumulvec(propotions = FALSE)` returns a vector of expected allocation counts, summing up to `n - exclude`. 
#' 

# For `pivec`,  an $M$-length vector with stationary/asymptotic visit frequencies.
#' 
#' 
#' @note When using the k-in-a-row design, set `matfun = kmatMarg` if using `pivec`, and otherwise `kmatFull`.
#'
#' @note At present, these functions are unable to incorporate in the calculations the impact of the recommended "fast start" stage for k-in-a-row and biased-coin designs. Such a stage begins with a classic UD run, until the first "minority" outcome is encountered (1 for below-median targets and vice versa). Generally such a fast start would make small-sample probability vectors approach the asymptotic distribution more quickly.
#' 
#' @inheritParams bcdmat
#' 
#' @param matfun The function to calculate the TPM. Depends on the specific design; see \code{\link{bcdmat}}. For all functions except `classicmat`, user must provide auxiliary parameters via `...`. For the k-in-a-row design, use `kmatMarg` for `pivec()` and `kmatFull` otherwise.
#' @param n For `currentvec, cumulvec`, at what step (= after how many observations) in the experiment would you like the vector calculated?
#' @param startdose (for `currentvec, cumulvec`), where does the experiment start? To be given as a dose-level index between 1 and \eqn{M}. If left as `NULL` (default), function will assume the equivalent of *"fair die roll"* among all doses. User can also specify your own \eqn{M}-length probability vector.
#' @param proportions Logical (`cumulvec` only) Would you like the results returned as proportions (= a probability vector; `TRUE`, default), or as cumulative allocation counts? 
#' @param exclude Integer (`cumulvec` only) Should the cumulative distribution exclude a certain number of initial observations? Default 0. 
#' @param marginalize logical (for `currentvec, cumulvec` when `matfun = kmatFull`) should the returned vector be marginalized over dose levels (`TRUE`, default), or should the full set with internal states be returned?
#' @param ... Arguments passed on to the design's matrix-calculating function.
#' 
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' @export
#'
#' @seealso
#'  - \code{\link{bcdmat}} for the functions calculating transition probability matrices for various up-and-down designs.
#'  - \code{\link{k2targ}} for target-finding design aids.
#'
#' 
#' @example inst/examples/vecExamples.r



#' @references 
#'  - Diaconis P, Stroock D. Geometric Bounds for Eigenvalues of Markov Chains. *Ann. Appl. Probab.* 1991;1(1):36-61. 
#'  - Hughes BD. *Random Walks and Random Environments, Vol. 1.* Oxford University Press, 1995.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50.

########################## Whew! Now the actual functions.


pivec <- function(cdf, matfun,  ...)
{
  checkCDF(cdf)
  if(identical(matfun, kmatFull)) stop("This function not set up to work with kmatFull(). Please kmatMarg().
  If you desire the full asymptotic distribution with internal states, use currentvec() with a large n.
                      5 * k * length(cdf) should suffice.")
  
    imat = matfun(cdf, ...) # matfun should take care of validations
  m = length(cdf)
  
# Pi is determined by the ratios between up/down movements,
#     represented as the TPM's two immediate off-diagonals
  vout = cumprod(c(1, imat[cbind(1:(m-1), 2:m)] / imat[cbind(2:m, 1:(m-1))] ))
# Normalizing to make it a probability vector
    return(vout / sum(vout))
}

######### Current allocation freq at patient n

#' @rdname pivec
#' @export
#' 
currentvec <- function(cdf, matfun, n, startdose = NULL, marginalize = TRUE, ...)
{
  requireNamespace('expm')
  if(identical(matfun, kmatMarg)) stop("For 'k-in-a-row' you must use kmatFull() rather than kmatMarg().\n")
  checkCDF(cdf)
  checkNatural(n, parname = 'n')
### Group UD! 
  if(identical(matfun, gudmat)) 
{
    message("For group UD, 'n' is interpreted as the number of cohorts rather than individual observations.\n")
}

# Grabbing the matrix now, in order to give...
mymat = matfun(cdf = cdf, ...)
m = nrow(mymat)

### Starting vector
# The NULL default sets a uniform starting vector

if(is.null(startdose)) { vec0=rep(1/m,m)
} else if (startdose %in% 1:m) {
    vec0=rep(0,m)
    vec0[startdose]=1
}
if(!exists('vec0')) vec0 = startdose
  # Now some validation 
if(any(vec0<0) || sum(vec0) != 1 || length(vec0)!=m) 
    stop("'startdose' must be a single dose level, or a probability vector over doses/states.\n")


### Then, the actual function could be just a one-liner :) 
#      We do need to de-matrix it though:
tmpout = as.numeric(vec0 %*% expm::`%^%` (mymat, n-1) )

# For non-krow designs, or if the full internal-state vector is desired for k-row - we're done.
if(!identical(matfun, kmatFull) | !marginalize) return(tmpout)

# For krow, we slog on to marginalize by dose: 
# First, need to figure out indirectly, whether low vs. high target
lowtarg = (mymat[1,2] > 0)
# Then k, while allowing for fluff or no fluff
kk = ceiling(m / length(cdf))
fluff = (m %% kk == 0)
cutts = seq(0, m+kk-1, kk) + 0.5
if(!lowtarg & !fluff) cutts = cutts - kk + 1

#return(tmpout)
sapply(split(tmpout, cut(1:m, cutts, labels = 1:length(cdf)) ), sum)
}

####### *Cumulative* allocation frequencies - perhaps the most practically useful

#' @rdname pivec
#' @export

cumulvec <- function(cdf, matfun, n, startdose = NULL, proportions = TRUE, exclude = 0, marginalize = TRUE, ...)
{
  requireNamespace('expm')
  if(identical(matfun, kmatMarg)) stop("For 'k-in-a-row' you must use kmatFull() rather than kmatMarg().\n")
  checkCDF(cdf)
  checkNatural(n, parname = 'n')
  checkNatural(exclude+1, parname = 'exclude', toolarge = n)
  ### Group UD! 
  if(identical(matfun, gudmat)) 
  {
    message("For group UD, 'n' is interpreted as the number of cohorts rather than individual observations.\n")
  }
  
  # Grabbing the matrix now, in order to give...
  progmat = matfun(cdf = cdf, ...)
  m = nrow(progmat)
  

### Starting vector
  # The NULL default sets a uniform starting vector
  if(is.null(startdose)) { vec0 = rep(1/m, m)
  } else if (startdose %in% 1:m) {
    vec0 = rep(0, m)
    vec0[startdose] = 1
  }
  if(!exists('vec0')) vec0 = startdose
  # Now some validation 
  if(any(vec0<0) || sum(vec0) != 1 || length(vec0) != m) 
    stop("'startdose' must be a single dose or a probability vector over doses.\n")
  
# The output vector  
  ovec = rep(0, m)
  for ( a in exclude:(n-1) ) ovec = ovec + vec0 %*% expm::`%^%` (progmat, a)
# Proportions or counts?
  if(proportions) ovec = ovec / (n - exclude)

# de-matrixing the output
  ovec = as.numeric(ovec)
  
  # For non-krow designs, or if the full internal-state vector is desired for k-row - we're done.
  if(!identical(matfun, kmatFull) | !marginalize) return(ovec)
  
  # For krow, we slog on to marginalize by dose: 
  # First, need to figure out indirectly, whether low vs. high target
  lowtarg = (progmat[1,2] > 0)
  # Then k, while allowing for fluff or no fluff
  kk = ceiling(m / length(cdf))
  fluff = (m %% kk == 0)
  cutts = seq(0, m+kk-1, kk) + 0.5
  if(!lowtarg & !fluff) cutts = cutts - kk + 1
  
  #return(tmpout)
  sapply(split(ovec, cut(1:m, cutts, labels = 1:length(cdf)) ), sum)
  

}





