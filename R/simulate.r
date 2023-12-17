############### Generalized Master Dose-Finding Function
# This function simulates sequential dose-finding experiments on a fixed dose grid			
# The response function is (implicitly) assumed monotone in [0,1]

# n: sample size
# target: target toxicity rate
# starting: the starting dose level
# cohort: the cohort size
# Fvals (vector or matrix): the true values of the response function on the dose grid
# progress: the function used to determine the next dose
# estimates: the function used to estimate/select the MTD
# progArgs: arguments needed for 'progress'. Need to be compatible for use in 'mapply'.
# thresholds (matrix): the response thresholds of participants, presented as percentiles rather than physical values. 
# nlev: the number of dose levels in the grid
# ensemble: the number of runs to be simulated


#' @export

dfsim <- function(n, starting=NULL, cohort=1, xvals, Fvals, progress = krow, progArgs = list(k=1), 
                  thresholds=NULL, seed = NULL,
              nlev=dim(Fvals)[1],ensemble=dim(Fvals)[2], quiet=FALSE, ...)
{
  require(plyr)

### Validation  
  if(nlev != length(xvals)) stop("Length mismatch (x values vs. F values).\n")
  if(var(sign(diff(xvals))) != 0) stop('Non-monotone x value set.\n')

  
  if(!quiet) cat(date(),'\n')		

  if(!is.null(seed)) set.seed(seed)
  est=NULL
  
  ## F values
  
  if (is.vector(Fvals)) {
    nlev=length(Fvals)
    Fvals=matrix(rep(Fvals,ensemble),nrow=nlev)
  }
  progArgs$maxlev=dim(Fvals)[1]
  
  #print(Fvals)
  if (is.null(thresholds)) thresholds=matrix(runif(n*ensemble),nrow=n)
  
  doses=matrix(NA,nrow=n+1,ncol=ensemble)
  responses=matrix(NA,nrow=n,ncol=ensemble)
  runid=1:ensemble
  
  if (is.null(starting)) {  	
    doses[1,]=sample(1:nlev,size=ensemble,replace=TRUE)
  } else doses[1,]=starting
  if (cohort>1) for (b in 2:cohort) doses[b,]=doses[1,]
  alive=1:ensemble
  
  for (a in seq(cohort+1,n+1,cohort))   ### main progression loop
  {
    ### We first obtain the current (a-1) patient's responses, then assign the next (a)
    ### Therefore, responses are only available up to patient a-1
    
    if(!quiet) cat('.')
    for (b in 1:cohort) responses[(a-b),alive]=ifelse(Fvals[cbind(doses[(a-b),alive],(1:ensemble)[alive])]>thresholds[(a-b),alive],1,0)
    
    alive=(!is.na(doses[a-1,]))
    #	cat(a-1,sum(alive),'\n')
    if(sum(alive)==0) break ### No more live runs; all have stopped
    #	cat(doses[a-1,alive],'\n')
  
    doses[a,alive]=mapply(FUN=progress,split(doses[1:(a-1),alive], col(matrix(doses[1:(a-1),alive], ncol=sum(alive))) ),
                split(responses[1:(a-1),alive], col(matrix(responses[1:(a-1),alive], ncol=sum(alive)))), MoreArgs=progArgs)
    # boundary conditions imposed by the Master
    doses[a,alive & doses[a,]>nlev]=nlev
    doses[a,alive & doses[a,]<1]=1
    
    ### in group designs, assigning the same to the remaining group:
    if (cohort>1 && a<=n) for (b in 2:min(cohort,n-a+1)) doses[a+b-1,]=doses[a,]
    if(a%%100==0) cat('\n')
    gc(verbose = FALSE)
  }
  
### Endgame

# "Dressing up" the dose levels (which are 1:m in the progress loop above) with real values
  doses = mapvalues(doses, 1:nlev, xvals)

  if(!quiet) cat('\n',date(),'\n')
  
  lout=list(scenarios=Fvals,sample=thresholds,dose=doses,response=responses,cohort=cohort,details=progArgs)
  return(lout)
}  ########  /dfsim

#------------------------- Implemented designs for dfsim()

#' @export

### k-in-a-row
krow <- function(doses, responses, k, hitarg=TRUE, cohort=1,fastStart=FALSE,...)
  ## firstpass: how many non-DLTs to require in the first pass? (i.e., the value of k for the first escalation) Ignored if 0.
{
  n=length(doses)
  dout=doses[n]
  
  if(hitarg)
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
bcd <- function(doses, responses, coin, hitarg=TRUE, fastStart=FALSE,...)
{
  n=length(doses)
  curr=doses[n]
  if(hitarg) {
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





