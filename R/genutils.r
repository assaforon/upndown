#' Up-and-Down Target Calculation and Calibration
#' 
#' Up-and-down target calculation and design options/guidance given a user-desired target.  
#' 
#' @details  
#' This suite of utilities helps users 
#' 
#'  - Figure out the approximate target response-rate, given design parameters
#'  - Suggest or specify design parameters, given user's target response-rate.
#'  
#'  Up-and-down designs (UDDs) generate random walks over dose space, with most dose-allocations usually taking place near the design's de-facto target percentile, called the **"balance point"** by some theorists to distinguish it from the user's designated target (Oron and Hoff 2009, Oron et al. 2022).
#'  
#'  Most k-in-a-row and group UDD parameter combinations yield balance points that are irrational percentiles of the dose-response function, and therefore are unappealing as official experimental targets.
#'  
#'  However, since the UD dose distribution has some width, and since even the balance point itself is only a close approximation for the actual average of allocated doses, the user's target **does not have to be identical to the balance point.** It only needs to be *"close enough"*.
#'  
#'  The `k2targ()` and `g2targ()` utilities are intended for users who already have a specific k-in-a-row or group design in mind, and only want to verify its balance point. The complementary utilities `ktargOptions(), gtargOptions()` provide a broader survey of design-parameter options within user-specified constraints, given a desired target.
#'  
#'  Lastly, `bcoin()` returns the biased-coin probabilities given the user's designated target. In contrast to the two other UDDs described above, the biased-coin design can target any percentile with a precisely matched balance point. That said, k-in-a-row and group UDDs offer some advantages over biased-coin in terms of properties and operational simplicity.


#------------------------------------- The actual functions ----------------------------#
### K-in-a-row functions

#' @export

k2targ<-function(k, hitarg=TRUE)
{
  checkNatural(k, parname = 'k', toolarge = 30)  
  tmp = 0.5 ^ (1/k)
  if(hitarg) return(tmp)
  1-tmp
}

#' @rdname k2targ
#' @export

ktargOptions<-function(target, tolerance = 0.1)
{
  if(length(target) > 1) stop("target must be a single number between 0 and 1.\n")
  checkTarget(target)
  
  hi = (target>=0.5) 
  if(!hi) target = 1-target
  
  klo = -1 / log2(target - tolerance)
  khi = -1 / log2(target + tolerance)
  krange = floor(klo):ceiling(khi)
  krange = krange[krange > 0]
  
  data.frame(k = krange, BalancePoint = k2targ(krange, hitarg = hi) )
}


## GUD functions

#' @rdname k2targ
#' @export


g2targ<-function(cohort, lower, upper)
{
checkNatural(c(cohort, lower+1, upper), parname = 'cohort, lower+1, upper', toolarge = 50)  
if(cohort<upper || upper<=lower) stop('Order must be lower < upper <= cohort.\n')

if(upper + lower == cohort) return(0.5)
  
uniroot(f=function(x, k, u, l) {pbinom(q=l, size=k, prob=x) + (pbinom(q=u-1, size=k, prob=x) - 1)}, interval=0:1, k=cohort, u=upper, l=lower)$root
}

# uniroot(f=function(x,kay,you,ell,bee1,bee2) {bee1*pbinom(q=ell,size=kay,prob=x)+bee2*(pbinom(q=you-1,size=kay,prob=x)-1)},interval=0:1,kay=k,you=u,ell=l,bee1=b1,bee2=b2)$root

#' @rdname k2targ
#' @export

gtargOptions<-function(target, maxsize = 6, tolerance = 0.1)
{
checkNatural(maxsize-1, parname = 'maxsize less 1', toolarge = 50)  
if(length(target) > 1) stop("target must be a single number between 0 and 1.\n")
checkTarget(target)

dout = data.frame()
for(k in 2:maxsize)
{
  refval = floor(target * k - 1e-10)
  for(l in 0:refval)
  {
    dtmp = data.frame()
    ucand = 2 * refval - l
    if(ucand > k) next
#    cat(k, l, ucand,'\n')
    u = ucand
   if(ucand > l) for(u in seq(ucand, l+1, -1))
    {
      bal = g2targ(cohort=k, lower=l, upper=u) 
      if(bal < target-tolerance) break
      if(bal > target+tolerance) next
      dtmp = rbind(dtmp, data.frame( cohort=k, lower=l, upper=u, BalancePoint=bal ) )
    }
    if(ucand < k) for(u in seq(max(l+1,ucand+1), k, 1))
    {
      bal = g2targ(cohort=k, lower=l, upper=u) 
      if(bal > target+tolerance) break
      if(bal < target-tolerance) next
      dtmp = rbind(dtmp, data.frame( cohort=k, lower=l, upper=u, BalancePoint=bal ) )
    }
    if(nrow(dtmp)>0) dout = rbind(dout, dtmp[order(dtmp$upper), ])
  }
}
return(dout)
}

################## BCD function

#' @rdname k2targ
#' @export

bcdCoin <- function(target, fraction = FALSE, nameplate = FALSE, tolerance = 0.02)
{
require(numbers)
  
## Validations
checkTarget(target)
checkTarget(tolerance, tname = "'tolerance'")
if(tolerance >= min(target, 1 - target)/2) stop("'tolerance' is set too large.\n")
if(tolerance < 1e-4) stop("'tolerance' is set unrealistically small. Change it to 1e-4 or more.\n")

### Targets too close to 0.5:

if(target >= 0.5 - tolerance && target <= 0.5 + tolerance) 
{
  cat("No need for coin with a target this close to 0.5. Just use Classical UD:
 - Move UP   after each negative response,                      
 - Move DOWN after each positive response.\n")
  
} else if(target < 0.5) { 
  
  coin = target / (1-target)
  cout = round(coin, digits = ceiling(-log10(tolerance)))
  if(fraction) 
# Find rational coin fraction, erring upwards b/c at <0.5 BCD biases a bit down
  {
    coinFrac = ratFarey(coin, n = round(1/tolerance), upper = TRUE) 
    if(nameplate)
    {
      frac2 = ratFarey(coin, n = round(1/tolerance), upper = FALSE) 
      if(abs(coin - frac2[1]/frac2[2]) < abs(coin - coinFrac[1]/coinFrac[2]))
        coinFrac = frac2
    }
    cout = paste(coinFrac[1], '/', coinFrac[2])
  }
  cat(paste("After positive response, move DOWN.
After negative response, 'toss a COIN':
   - with probability of approximately", cout, 'move UP
   - Otherwise REPEAT same dose.\n') )
  
} else if(target > 0.5) { 
  
  coin = (1-target) / target 
  cout = round(coin, digits = ceiling(-log10(tolerance)))
  if(fraction) 
    # Find rational coin fraction, erring upwards b/c at <0.5 BCD biases a bit down
  {
    coinFrac = ratFarey(coin, n = round(1/tolerance), upper = FALSE) 
    if(nameplate)
    {
      frac2 = ratFarey(coin, n = round(1/tolerance), upper = TRUE) 
      if(abs(coin - frac2[1]/frac2[2]) < abs(coin - coinFrac[1]/coinFrac[2]))
        coinFrac = frac2
    }
    cout = paste(coinFrac[1], '/', coinFrac[2])
  }
  cat(paste("After negative response, move UP.
After positive response, 'toss a COIN':
   - with probability of approximately", cout, 'move DOWN
   - Otherwise REPEAT same dose.\n') )
}  
  
}



############################### Auxiliary validation utilities


validUDinput<-function(cdf,target)
{
  checkCDF(cdf)
  checkTarget(target)
  if(length(cdf) < 3) stop ("These designs don't work with <3 dose levels.\n")
}

checkTarget <- function(target, tname = 'Target')
  if(target<=0 || target>=1) stop(paste(tname, "has to be in (0,1).\n"))

checkCDF <- function(cdf)
  if(min(cdf)<0 || max(cdf)>1 || any(diff(cdf) < 0) || var(cdf)==0) stop("cdf should be a CDF.\n")

## Natural number verification

checkNatural <- function(k, parname, toolarge = 1000)
if (any(k != round(k) | k < 1 | k >= toolarge)) 
  stop(parname, "must be a natural number smaller than ", toolarge, ".\n")

