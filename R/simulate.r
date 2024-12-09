##' Generalized Dose-Finding Ensemble Simulator
#'
#' This function simulates sequential dose-finding experiments on a fixed dose grid. The response function is (implicitly) assumed monotone in `(0,1)`
#' 
#' A vectorized dose-finding simulator, set up to run an entire ensemble simultaneously. 
#' The simulated doses are indices `1:nlev` with `nlev` being the number of dose levels.
#' Upon output they can be optionally "dressed up" with physical values using the `xvals` argument.
#' 
#' The simulator's essential use within the `upndown` package is to estimate bootstrap confidence intervals for dose-averaging target estimates, via \code{\link{dfboot}}. But it can be also used stand-alone as a study-design aid.
#' 
#' The particular dose-finding design simulated is determined by `design` and its argument list `desArgs`. The 3 straightforward extensions of the median-finding "Classical" UDD are available, namely "k-in-a-row", biased-coin and group (cohort) UDD. To simulate the the median-finding "Classical" UDD itself, use `krow` with `desArgs = list(k=1)`. 
#' Other non-UDD dose-finding designs - e.g., CRM, CCD, BOIN, etc. - can also be made compatible with `dfsim`. Utilities to run those 3 in particular are available on GitHub, under `assaforon/UpndownBook/P3_Practical/OtherDesigns.r`.
#' 
#' If you want to create a `design` function yourself, it would need to accept `doses, responses` as input, and return the next dose allocation (as an integer index).
#' The main progression loop is run via `mapply`.

#' @param n sample size
#' @param starting the starting dose level. If `NULL` (default), will be randomized.
#' @param sprobs the probability weights if using a randomized starting dose. If `NULL` (default) will be discrete-uniform.
#' @param cohort the cohort (group) size, default 1.
#' @param Fvals (vector or matrix): the true values of the response function on the dose grid. These are the dose-response scenarios from which the experimental runs will be simulated. If running an ensemble with different scenarios, each scenarios is a column. If running an identical-scenario ensemble, provide a single vector as well as `ensemble`.
#' @param ensemble the number of different runs/scenarios to be simulated. Will be determined automatically if `Fvals` is a matrix, as the number of columns.
#' @param design the dose-finding design function used to determine the next dose. Default `krow`; see \code{\link{krow}} for options.
#' @param desArgs List of arguments passed on to `design`. Need to be compatible for use in `mapply`. Default is `list(k=1)`, which together with `design = krow` will generate a Clasical (median-finding) UDD simulation.
#' @param thresholds Matrix of size (at least) `n` by `ensemble`, the response thresholds of participants, presented as percentiles (i.e., output of `runif()`) rather than physical values. If `NULL` (default), they will be simulated on the fly. When running comparative performance simulations, we recommend providing the same thresholds to everything you want to compare on equal footing.
#' @param seed The random seed if simulating the thresholds. Can be kept *"floating"* (i.e., varying between calls) if left as `NULL` (default).
#' @param showdots Logical: print out a dot (`.`) after each designion step in `1:n`, and the start/end time stamps? Default `TRUE`.
#' 
#' @author Assaf P. Oron
#' 
#' @note This is an adaptation of a non-package function used by the author for well over a decade before incorporating it into `upndown` in late 2023. For initial guidance, see the code example. If you encounter any funny behavior please let me know. Thank you!

#' @seealso 
#' - \code{\link{dfboot}}

#' @example inst/examples/simExample.r


#' @return A list with the following elements:
#'  - `scenarios`: `Fvals`
#'  - `sample`: `thresholds`
#'  - `doses`: The matrix of simulated dose allocations for each run (`n+1` by `ensemble`)
#'  - `responses`: The matrix of simulated responses (0 or 1) for each run (`n` by `ensemble`)
#'  - `cohort`: `cohort`
#'  - `details`: `desArgs`


#' @export

dfsim <- function(n, starting=NULL, sprobs = NULL, cohort=1, Fvals, ensemble = dim(Fvals)[2], 
                  design = krow, desArgs = list(k=1), thresholds=NULL, seed = NULL, showdots = TRUE)
{

### Validation  

 checkNatural(n)
 ## F values
  
  if (is.vector(Fvals)) {
    checkCDF(Fvals)
    nlev=length(Fvals)
    Fvals=matrix(rep(Fvals,ensemble),nrow=nlev)
  } else if(length(dim(Fvals)) == 2) {
    apply(Fvals, 2, checkCDF)
    nlev = dim(Fvals)[1]
  } else stop("'Fvals' must be a vector or matrix.\n")
  
 desArgs$maxlev = nlev
 
  if(showdots) cat(date(),'\n')		
  if(!is.null(seed)) set.seed(seed)

  #print(Fvals)

###### Prep

## Response threholds 
   if (is.null(thresholds)) thresholds=matrix(runif(n*ensemble),nrow=n) else {
    
    checkTarget(as.vector(thresholds))
    if(nrow(thresholds) < n || ncol(thresholds) < ensemble) stop("Not enough random thresholds.\n")
  }
  
  doses=matrix(NA,nrow=n+1,ncol=ensemble)
  responses=matrix(NA,nrow=n,ncol=ensemble)
  runid=1:ensemble
  
## Randomized starting dose (in case startdose not given)
  if (is.null(starting)) {  
    if(is.null(sprobs)) sprobs = rep(1/nlev, nlev)
    doses[1,]=sample(1:nlev,size=ensemble,replace=TRUE, prob = sprobs)
  } else doses[1,]=starting
  if (cohort>1) for (b in 2:cohort) doses[b,]=doses[1,]
  alive=1:ensemble
  
###-------------------------- main progression loop -----------------------###

  for (a in seq(cohort+1,n+1,cohort))  
  {
    ### We first obtain the current (number a-1) responses, then assign the next (a)
    ### Therefore, responses are only available up to a-1
    
    if(showdots) cat('.')
    for (b in 1:cohort) responses[(a-b),alive]=ifelse(Fvals[cbind(doses[(a-b),alive],(1:ensemble)[alive])]>thresholds[(a-b),alive],1,0)
    
    alive=(!is.na(doses[a-1,]))
    #	cat(a-1,sum(alive),'\n')
    if(sum(alive)==0) break ### No more live runs; all have stopped
    #	cat(doses[a-1,alive],'\n')
  
    doses[a,alive]=mapply(FUN=design, split(doses[1:(a-1),alive], col(matrix(doses[1:(a-1),alive], ncol=sum(alive))) ),
                split(responses[1:(a-1),alive], col(matrix(responses[1:(a-1),alive], ncol=sum(alive)))), MoreArgs=desArgs)
    # boundary conditions imposed by the Master
    doses[a,alive & doses[a,]>nlev]=nlev
    doses[a,alive & doses[a,]<1]=1
    
    ### in group designs, assigning the same to the remaining group:
    if (cohort>1 && a<=n) for (b in 2:min(cohort,n-a+1)) doses[a+b-1,]=doses[a,]
    if(a%%100==0) cat('\n')
    gc(verbose = FALSE)
  }
  
####### Endgame

  if(showdots) cat('\n',date(),'\n')
  
  lout=list(scenarios=Fvals, sample=thresholds, doses=doses, responses=responses, cohort=cohort, details=desArgs)
  return(lout)
}  ########  /dfsim


#------------------------- In-package implemented designs for dfsim()

#' Up-and-Down Design Rules for use in Dose-Finding Simulator
#' 
#' Rules for k-in-a-row, Biased-Coin UD, and Group UD, coded as functions compatible with
#'    the generic dose-finding simulator `dfsim()`
#'    
#' These functions work on each virtual experimental run individually. 
#' 
#' 
#' Rules for some popular or well-studied non-up-and-down
#'    
#'  
#' @inheritParams bcdmat
#'   
#' @param doses,responses (mandatory arguments) vectors of the run's current sequence of doses (in ordinal/index scale) and responses
#' @param lowTarget (`krow` and `bcd`) logical: is the target below 0.5 (median threshold)? 
#' @param fastStart (`krow` and `bcd`) logical: should the experiment begin with a classical-UD-like stage until the first "minority" response is observed (i.e., a 1 for below-median targets and vice versa)? Even though `TRUE` delivers better experimental performance and is recommended when allowed, default is `FALSE` because toxicity/safety studies are unlikely to allow it. 
#' @param coin (`bcd` only) the biased-coin probability. Note that unlike `bcdmat()`, here the function does not figure out automatically the coin probability and upper/lower target location from the provided target. 
#' @param s (`groupUD` only) the group/cohort size, analogous to `cohort` in `gudmat()`. We use a different name here because `cohort` is already used in `dfsim` that calls these utilities.
#' @param ll,ul (`groupUD` only) how many positive responses are allowed for a move upward, and how many are required for a move downward, respectively. Analogous to `lower, upper` in `gudmat()`. For example `s=3, ll=0, ul=2` evaluates groups of 3 observations at a time, moves up if none are positive, down if \eqn{>=2} are positive, and repeats the same dose with 1 positive.
#' @param ... Technical pass-through argument, to allow for flexibility when constructing design-comparison simulation ensembles.

#' @return the next dose allocation

#' @export

### k-in-a-row
krow <- function(doses, responses, k, lowTarget=NULL, fastStart=FALSE,...)
{
  if(is.null(lowTarget)) if(k>1) stop('Must provide `lowTarget`!\n') else lowTarget = FALSE
  n=length(doses)
  dout=doses[n]
  
  if(!lowTarget)
  {
    if(responses[n]==0) return(dout+1) 
    if(fastStart && sum(responses)==n) return(dout-1)
    if(n<k) return(dout)
    if(mean(responses[(n-k+1):n])==1 && (k==1 || var(doses[(n-k+1):n])==0)) return(dout-1)
    return(dout)
  } else {
    
    if(fastStart && sum(responses)==0) return(dout+1)
    if(responses[n]==1) return(dout-1)  # toxicity in current observation: down
    if(n<k) return(dout)
    if(sum(responses[(n-k+1):n])==0 && (k==1 || var(doses[(n-k+1):n])==0)) return(dout+1)
    return(dout)
  }
}

#' @inheritParams bcdmat
#' @rdname krow
#' @export

### BCD
bcd <- function(doses, responses, coin, lowTarget, fastStart=FALSE,...)
{
  n=length(doses)
  curr=doses[n]
  if(!lowTarget) {
    dout=ifelse(responses[n]==0,curr+1,ifelse(runif(1)<=coin,curr-1,curr))
    if(fastStart && sum(responses)==n) dout=curr-1
  } else {
    dout=ifelse(responses[n]==1,curr-1,ifelse(runif(1)<=coin,curr+1,curr))
    if(fastStart && sum(responses)==0) dout=curr+1
  }
  return(dout)
}

#' @rdname krow
#' @export

### Group UD
groupUD=function(doses, responses, s, ll, ul,...)
{
  if(ll>s || ul>s) stop('Group up-down boundaries cannot be greater than group size.\n')
  if (ll>=ul) stop('ll bound cannot be greater than ul bound.\n')
  
  n=length(doses)
  curr=doses[n]
  if(n%%s>0) return(curr) # only evaluating when group is full
  if(n<s) return(curr)
  dlt=sum(responses[(n-s+1):n])
  dout=ifelse(dlt<=ll, curr+1, ifelse(dlt>=ul, curr-1, curr))
  return(dout)
}





