#' Up-and-Down Target Calculation and Design Guidance
#' 
#' Up-and-down target calculation, as well as design options/guidance given a user-desired target.  
#' 
#' @details  
#' This suite of utilities helps users 
#' 
#'  - Figure out the approximate target response-rate (a.k.a. the *balance point*), given design parameters;
#'  - Suggest potential design parameters, given user's desired target response-rate and other constraints.
#'  
#'  Up-and-down designs (UDDs) generate random walks over dose space, with most dose-allocations usually taking place near the design's de-facto target percentile, called the **"balance point"** by some theorists to distinguish it from the user-designated target in case they differ (Oron and Hoff 2009, Oron et al. 2022).
#'  
#'  Most k-in-a-row and group UDD parameter combinations yield balance points that are irrational percentiles of the dose-response function, and therefore are unappealing as official experimental targets.
#'  
#'  However, since the UD dose distribution has some width, and since even the balance point itself is only a close approximation for the actual average of allocated doses, the user's target **does not have to be identical to the balance point.** It only needs to be *"close enough"*.
#'  
#'  The `k2targ()` and `g2targ()` utilities are intended for users who already have a specific k-in-a-row or group design in mind, and only want to verify its balance point. The complementary utilities `ktargOptions(), gtargOptions()` provide a broader survey of design-parameter options within user-specified constraints, given a desired target. The 0.05 default tolerance level around the target, is what we recommend as *"close enough"*. Otherwise, it is probably better to use biased-coin.
#'  
#'  Lastly, `bcoin()` returns the biased-coin probabilities given the user's designated target. In contrast to the two other UDDs described above, the biased-coin design can target any percentile with a precisely matched balance point. That said, k-in-a-row and group UDDs offer some advantages over biased-coin in terms of performance and operational simplicity. `bcoin()` can return the probability as a decimal, rational fraction, or both (default). 
#'  
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	
#'  
#' @seealso
#'  - \code{\link{bcdmat}} for the functions calculating transition probability matrices for various up-and-down designs.
#'  - \code{\link{pivec}} for functions calculating key probability vectors for the designs.
#'
#' 
#' @param target the desired target response rate (as a fraction in \eqn{(0,1)}), where relevant.
#' @param k the number of consecutive identical responses required for dose transitions (k-in-a-row functions only).
#' @param lowTarget logical, `k2targ()` only: is the design targeting below-median percentiles, with \eqn{k} repeated negative responses needed to level up and only one to level down - or vice versa? 
#' @param cohort,lower,upper `g2targ()` only: the cohort (group) size, how many positive responses are allowed for a move upward, and how many are required for a move downward, respectively. For example `cohort=3, lower=0, upper=2` evaluates groups of 3 observations at a time, moves up if none are positive, down if \eqn{>=2} are positive, and repeats the same dose with 1 positive.
#' @param tolerance 
#'  - For `ktargOptions(), gtargOptions()`: the half-width of the interval around `target` in which to search for design options. Default 0.05. 
#'  - For `bcoin()`: the half-width of the interval around 0.5 in which the function recommends to simply use classical UD without a coin. Default 0.05.
#' @param maxk `ktargOptions()` only: the maximum value of \eqn{k} to consider.
#' @param minsize,maxsize `gtargOptions()` only: the minimum and maximum cohort size to consider. `minsize` has to be at least 2 (cohort size 1 is equivalent to classical UD).
#' @param presentation `bcoin()` only: whether to report the coin probability as a `'decimal'`,  rational `'fraction'`. or `'both'` (default). 
#' @param digits how many digits to round to? Default 4. Since the functions use R's `round()` utility, trailing zeroes will be truncated.
#' 
#' 
#' @return 
#'  - `k2targ(), g2targ()`: the official balance point given the user-provided design parameters.
#'  - `ktargOptions(), gtargOptions()`: a `data.frame` with design parameters and balance points, for all options that meet user-provided constraints. A printed string provides dose transition rule guidance.
#'  - `bcoin():` a printed string that informs user of the biased-coin design rules, including the 'coin' probability in its user-chosen format (decimal or fraction). In case the user-desired target is 0.5 or very close to it, the string will inform user that they are better off just using classical UDD without a coin.
#'  
#'  
#' @references 
#'  - Durham SD, Flournoy N. Random walks for quantile estimation. In: *Statistical Decision Theory and Related Topics V* (West Lafayette, IN, 1992). Springer; 1994:467-476.
#'  - Gezmu M, Flournoy N. Group up-and-down designs for dose-finding. *J Stat Plan Inference.* 2006;136(6):1749-1764.
#'  - Oron AP, Hoff PD. The k-in-a-row up-and-down design, revisited. *Stat Med.* 2009;28:1805-1820.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50.

#------------------------------------- The actual functions ----------------------------#
###-------------- K-in-a-row functions

#' @export

k2targ<-function(k, lowTarget, digits = 4)
{
  checkNatural(k, parname = 'k', toolarge = 30)  
  tmp = 0.5 ^ (1/k)
  if(lowTarget) return(round(1-tmp, digits) )
  round(tmp, digits)
}

#' @rdname k2targ
#' @export

ktargOptions<-function(target, tolerance = 0.05, maxk = 20, digits = 4)
{
  if(length(target) > 1) stop("target must be a single number between 0 and 1.\n")
  checkTarget(target)
  checkTarget(target - tolerance, tname = "target minus tolerance")

# For simplicity, doing the above-median arithmetic
  hi = (target>=0.5) 
  if(!hi) target = 1-target
  
  hiend = k2targ(maxk, lowTarget = FALSE)
  
  klo = -1 / log2( max(0.5, target - tolerance) )
  khi = -1 / log2( min(hiend, target + tolerance) )
  
#  return(c(klo, khi))
  krange = ceiling(klo):floor(khi)
  if(any(diff(krange) < 0 ) ) return("No k options within target range; please increase tolerance band or use biased-coin.")
  krange = sort( krange[krange > 0 & krange <= maxk] )
  
  message("For the following targets, k consecutive ", ifelse(hi, "positive", "negative"),
      " responses are needed for a move ", ifelse(hi, "DOWN.\n", "UP.\n"),
      "Only one ", ifelse(hi, "negative", "positive"), " response is needed for the opposite move.\n")
  
  data.frame(k = krange, BalancePoint = k2targ(krange, lowTarget = !hi, digits = digits) )
}


##--------------------- GUD functions

#' @rdname k2targ
#' @import stats
#' @export


g2targ<-function(cohort, lower, upper, digits = 4)
{
checkNatural(c(cohort, lower+1, upper), parname = 'cohort, lower+1, upper', toolarge = 50)  
if(cohort<upper || upper<=lower) stop('Order must be lower < upper <= cohort.\n')

if(upper + lower == cohort) return(0.5)
  
round(uniroot(f=function(x, k, u, l) {pbinom(q=l, size=k, prob=x) + (pbinom(q=u-1, size=k, prob=x) - 1)}, interval=0:1, k=cohort, u=upper, l=lower)$root, digits)
}

# uniroot(f=function(x,kay,you,ell,bee1,bee2) {bee1*pbinom(q=ell,size=kay,prob=x)+bee2*(pbinom(q=you-1,size=kay,prob=x)-1)},interval=0:1,kay=k,you=u,ell=l,bee1=b1,bee2=b2)$root

#' @rdname k2targ
#' @export

# Vals
gtargOptions<-function(target, minsize = 2, maxsize = 6, tolerance = 0.05, digits = 4)
{
if(maxsize<minsize) stop("'maxsize' cannot be smaller than 'minsize'.\n")
checkNatural(c(maxsize-1, minsize-1), parname = 'Cohort size less 1', toolarge = 50)  
if(length(target) > 1) stop("target must be a single number between 0 and 1.\n")
checkTarget(target)
# \vals

dout = data.frame()
for(k in minsize:maxsize)
{
  refval = floor(target * k - 1e-10)
  for(l in 0:refval)
  {
    dtmp = data.frame()
    ucand = 2 * refval - l
    if(ucand > k) next
#    message(k, l, ucand,'\n')
    u = ucand
   if(ucand > l) for(u in seq(ucand, l+1, -1))
    {
      bal = g2targ(cohort=k, lower=l, upper=u, digits = digits) 
      if(bal < target-tolerance) break
      if(bal > target+tolerance) next
      dtmp = rbind(dtmp, data.frame( Cohort=k, Lower=l, Upper=u, BalancePoint=bal ) )
    }
    if(ucand < k) for(u in seq(max(l+1,ucand+1), k, 1))
    {
      bal = g2targ(cohort=k, lower=l, upper=u, digits = digits) 
      if(bal > target+tolerance) break
      if(bal < target-tolerance) next
      dtmp = rbind(dtmp, data.frame( Cohort=k, Lower=l, Upper=u, BalancePoint=bal ) )
    }
    if(nrow(dtmp)>0) dout = rbind(dout, dtmp[order(dtmp$Upper), ])
  }
}

if(nrow(dout) == 0) return("No Group UDD options within target range; please increase tolerance band or use biased-coin.")

message("For each design, if positive responses <= Lower, move UP.\n")
message("                 if positive responses >= Upper, move DOWN.\n")
message("  otherwise, REPEAT same dose (relevant only when Upper - Lower > 1).\n")

return(dout)
}

################## BCD function

#' @rdname k2targ
#' @export

bcoin <- function(target, presentation = 'both', digits = 4, tolerance = 0.05)
{
requireNamespace('MASS')
  
## Validations
checkTarget(target)
checkTarget(tolerance, tname = "'tolerance'")
if(tolerance >= min(target, 1 - target)/2) stop("'tolerance' is set too large.\n")
if(tolerance < 1e-4) stop("'tolerance' is set unrealistically small. Change it to 1e-4 or more.\n")
allowedprez = c('decimal', 'fraction', 'both') 
if(!(presentation %in% allowedprez) ) 
  stop("Presentation must be one of ", paste(allowedprez, collapse = ', ') )

### Targets too close to 0.5:

if(target >= 0.5 - tolerance && target <= 0.5 + tolerance) 
{
  message("No need for coin with a target this close to 0.5. Just use Classical UD:
 - Move UP   after each negative response,                      
 - Move DOWN after each positive response.\n")
  
} else if(target < 0.5) { 
  
  coin = target / (1-target)
  if(presentation == 'decimal') cout = round(coin, digits) else
    {
      frac = MASS::fractions(coin, cycles = 4)
      cout = frac
      if(presentation == 'both') cout = paste(cout, '  (', round(coin, digits), ')', sep='' )
    }
  
  message(paste("After a positive response, move DOWN.
After a negative response, 'toss a Coin':
   - with probability of", cout, 'move UP
   - Otherwise REPEAT same dose.\n') )
  
} else if(target > 0.5) { 
  
  coin = (1-target) / target 
  if(presentation == 'decimal') cout = round(coin, digits) else
  {
    frac = MASS::fractions(coin, cycles = 4)
    cout = frac
    if(presentation == 'both') cout = paste(cout, '  (', round(coin, digits), ')', sep='' )
  }
  
  message(paste("After a negative response, move UP.
After a positive response, 'toss a Coin':
   - with probability of", cout, 'move DOWN
   - Otherwise REPEAT same dose.\n') )
}  
  
}



############################### Auxiliary validation utilities

#' Data Validation Utilities for `upndown`
#' 
#' Validation of input values
#' 
#' 
#' @param cdf vector of values, should be nondecreasing between 0 and 1 (inclusive)
#' @param target numeric value(s), should be between 0 and 1 (exclusive)
#' @param parname,tname string, name of variable to plug in for reporting the error back
#' @param k (`checkNatural()` only) input number to check whether it's a natural number 
#' @param toolarge (`checkNatural()` only) what number would be considered too large to be realistic?
#' @param x (`checkDose()` only) input object to be verified as valid dose values
#' @param maxfrac (`checkDose()` only) maximum number of unique values (as fraction of sample size) considered realistic for up-and-down data. Default \eqn{0.9n}. Function also gives a warning if number exceeds \eqn{n/2}, since typically this suggests the effective sample size around target is too small.
#' @param y (`checkResponse()` only) input object to be verified as valid response values (`TRUE/FALSE or 0/1)
#' @param flatOK logical (`checkCDF()` only) if the CDF is completely flat, should function issue a warning (`TRUE`, default) - or stop with an error (`FALSE`)?
#' 
#' 
#' @return If a validation issue is found, these functions stop with a relevant error message. If no issue is found, they run through without returning a value.

#' @export
validUDinput<-function(cdf, target)
{
  checkCDF(cdf)
  checkTarget(target)
  if(length(cdf) < 3) stop ("These designs don't work with <3 dose levels.\n")
}

#' @rdname validUDinput
#' @export
checkTarget <- function(target, tname = 'Target')
  if( any(target<=0 | target>=1) ) stop(paste(tname, "has to be in (0,1).\n"))

#' @rdname validUDinput
#' @import stats
#' @export
checkCDF <- function(cdf, flatOK = TRUE)
{
  if(length(cdf) < 1 || min(cdf)<0 || max(cdf)>1 || any(diff(cdf) < 0)) stop("cdf should be a CDF.\n")
  if(var(cdf)==0) {
    if(!flatOK) stop("cdf is compeletely flat.\n")
    warning("cdf is compeletely flat.\n")
  }
}

## Natural number verification

#' @rdname validUDinput
#' @export
checkNatural <- function(k, parname, toolarge = 1000)
if (any(k != round(k) | k < 1 | k >= toolarge)) 
  stop(parname, "must be a natural number smaller than ", toolarge, ".\n")

## Response coding

#' @rdname validUDinput
#' @export
checkDose <- function(x, maxfrac = 0.9)
{
  if(length(dim(x))>1 || any(dim(x))>1) stop("Dose must be a vector or equivalent.\n")
  if( length(unique(x)) > maxfrac * length(x))
    stop("Too many unique dose values. Experiment not U&D, or was too short, or data-quality error.\n")
  if( length(unique(x)) > length(x)/2)
    warning("Number of unique dose values exceeds n/2.\n The experiment might be too short, or started too far from apparent target.\n Otherwise, there might be a data-quality error.\n")
}

#' @rdname validUDinput
#' @export
checkResponse <- function(y)
{
if(length(dim(y))>1 || any(dim(y))>1) stop("Response must be a vector or equivalent.\n")
if( any( !is.logical(y) & !(y %in% 0:1) ) )
  stop("Response must be logical or coded as 0/1.\n")
}

