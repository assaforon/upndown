##' Generalized Dose-Finding Ensemble Simulator
#'
#' This function simulates sequential dose-finding experiments on a fixed dose grid. The response function is (implicitly) assumed monotone in `(0,1)`
#' 
#' A vectorized dose-finding simulator, set up to run an entire ensemble simultaneously. 
#' The simulated doses are indices `1:nlev` with `nlev` being the number of dose levels.
#' Upon output they can be optionally "dressed up" with physical values using the `xvals` argument.
#' 
#' The simulator's essential use within the `upndown` package is to estimate confidence intervals for dose-averaging target estimates. But it can be also used stand-alone as a study-design aid.
#' 
#' The particular dose-finding design simulated is determined by `progress` and its argument list `progArgs`. UDD design functions are provided, but other designs - CRM, CCD, etc. - are also compatible. CRM and CCD in particular are available from the author upon request.
#' The `progress` functions need to accept `doses, responses` as input, and return the next dose allocation (as an index).
#' The progression loop is run via `mapply`.

#' @param n sample size
#' @param starting the starting dose level. If `NULL` (default), will be randomized.
#' @param sprobs the probability weights if using a randomized starting dose. If `NULL` (default) will be discrete-uniform.
#' @param cohort the cohort (group) size, default 1.
#' @param Fvals (vector or matrix): the true values of the response function on the dose grid. These are the dose-response scenarios from which the experimental runs will be simulated. If running an ensemble with different scenarios, each scenarios is a column. If running an identical-scenario ensemble, provide a single vector as well as `nlev, ensemble`.
#' @param nlev The number of dose levels in the grid. Will be determined automatically if `Fvals` is a matrix, as the number of rows. 
#' @param ensemble the number of different runs/scenarios to be simulated. Will be determined automatically if `Fvals` is a matrix, as the number of columns.
#' @param progress the dose-finding design function used to determine the next dose. Default `krow`; see \code{\link{krow}} for options.
#' @param progArgs List of arguments passed on to `progress`. Need to be compatible for use in `mapply`. Default is `list(k=1)`, which together with `progress = krow` will generate a Clasical (median-finding) UDD simulation.
#' @param thresholds Matrix of size (at least) `n` by `ensemble`, the response thresholds of participants, presented as percentiles (i.e., output of `runif()`) rather than physical values. If `NULL` (default), they will be simulated on the fly
#' @param seed The random seed if simulating the thresholds. Can be kept *"floating"* (i.e., varying between calls) if left as `NULL` (default).
#' @param quiet Logical: suppress printing out a dot (`.`) after each progression step in `1:n`, and the start/end time stamps? Default `FALSE`.
#' 
#' @author Assaf P. Oron
#'
#' @note This is an adaptation of a non-package function used by the author for well over a decade before incorporating it into `upndown` in late 2023. If you encounter any funny behavior please let me know. Thank you!


#' @return A list with the following elements:
#'  - `scenarios`: `Fvals`
#'  - `sample`: `thresholds`
#'  - `dose`: The matrix of simulated dose allocations for each run (`n+1` by `ensemble`)
#'  - `response`: The matrix of simulated responses (0 or 1) for each run (`n` by `ensemble`)
#'  - `cohort`: `cohort`
#'  - `details`: `progArgs`

#' @export

dfsim <- function(n, starting=NULL, sprobs = NULL, cohort=1, Fvals, nlev=dim(Fvals)[1], ensemble = dim(Fvals)[2], 
                  progress = krow, progArgs = list(k=1), thresholds=NULL, seed = NULL, quiet = FALSE)
{

### Validation  

 checkNatural(n)
 ## F values
  
  if (is.vector(Fvals)) {
    checkCDF(Fvals)
    nlev=length(Fvals)
    Fvals=matrix(rep(Fvals,ensemble),nrow=nlev)
  } else apply(Fvals, 2, checkCDF)
 progArgs$maxlev=dim(Fvals)[1]
 
  if(!quiet) cat(date(),'\n')		
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
    
    if(!quiet) cat('.')
    for (b in 1:cohort) responses[(a-b),alive]=ifelse(Fvals[cbind(doses[(a-b),alive],(1:ensemble)[alive])]>thresholds[(a-b),alive],1,0)
    
    alive=(!is.na(doses[a-1,]))
    #	cat(a-1,sum(alive),'\n')
    if(sum(alive)==0) break ### No more live runs; all have stopped
    #	cat(doses[a-1,alive],'\n')
  
    doses[a,alive]=mapply(FUN=progress, split(doses[1:(a-1),alive], col(matrix(doses[1:(a-1),alive], ncol=sum(alive))) ),
                split(responses[1:(a-1),alive], col(matrix(responses[1:(a-1),alive], ncol=sum(alive)))), MoreArgs=progArgs)
    # boundary conditions imposed by the Master
    doses[a,alive & doses[a,]>nlev]=nlev
    doses[a,alive & doses[a,]<1]=1
    
    ### in group designs, assigning the same to the remaining group:
    if (cohort>1 && a<=n) for (b in 2:min(cohort,n-a+1)) doses[a+b-1,]=doses[a,]
    if(a%%100==0) cat('\n')
    gc(verbose = FALSE)
  }
  
####### Endgame

  if(!quiet) cat('\n',date(),'\n')
  
  lout=list(scenarios=Fvals,sample=thresholds,dose=doses,response=responses,cohort=cohort,details=progArgs)
  return(lout)
}  ########  /dfsim


#------------------------- Implemented designs for dfsim()

#' @export

### k-in-a-row
krow <- function(doses, responses, k, lowTarget=NULL, cohort=1, fastStart=FALSE,...)
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
    if(sum(responses[(n-cohort+1):n])>=1) return(dout-1)  # toxicity in current observation/cohort: down
    if(n<k) return(dout)
    if(sum(responses[(n-k+1):n])==0 && (k==1 || var(doses[(n-k+1):n])==0)) return(dout+1)
    return(dout)
  }
}

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

#' @export

### Group UD
groupUD=function(doses,responses,s,ll,ul,...)
{
  if(ll>s || ul>s) stop('Group up-down boundaries cannot be greater than group size.\n')
  if (ll>=ul) stop('Lower bound cannot be greater than upper bound.\n')
  
  n=length(doses)
  curr=doses[n]
  if(n%%s>0) return(curr) # only evaluating when group is full
  if(n<s) return(curr)
  dlt=sum(responses[(n-s+1):n])
  dout=ifelse(dlt<=ll, curr+1, ifelse(dlt>=ul, curr-1, curr))
  return(dout)
}





