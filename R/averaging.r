### Up-and-Down target estimates based on dose averaging 


##' Original Dixon and Mood (1948) point estimate
#'
##' Basic version; formula assumes uniform spacing but should work anyway
#'
#' In their documentation of the Up-and-Down algorithm, Dixon and Mood (1948) presented an estimation method based on tallying responses

#' @inheritParams reversmean
#' 
#' @param flip logical: should we flip D-M's approach and use the more-common outcome? (default \code{FALSE})

dixonmood<-function(x, y, full=FALSE, flip=FALSE)
{
n=length(x)
if (!(length(y)==n)) stop('Mismatched lengths.\n')
if(length(unique(y))==1) return(NA) # degenerate case, all 0's or 1's
xvals=sort(unique(x))
spacing=mean(diff(xvals))
# D-M method uses the less frequent outcome, hence might drop leading/trailing "tails"
chosen=ifelse(mean(y)<0.5,1,0)
if(flip) chosen=1-chosen
track=x[y==chosen]
# After all this prep, the formula is anticlimactic :)
if(full) return(c(A=sum((track-min(track))/spacing),N=length(track),d=spacing))
mean(track)+spacing*(0.5-chosen)
}

#----------------------- Utilities for reversal-anchored averaging -----------------#




#' Reversal-anchored averaging estimators for Up-and-Down
#' 
#' Dose-averaging target estimation for Up-and-Down experiments, historically the most popular approach, but not recommended as primary nowadays. Provided for completeness.
#' 
#' Up-and-Down designs (UDDs) allocate doses in a random walk centered nearly symmetrically around a balance point. Therefore, a modified average of allocated doses could be a plausible estimate of the balance point's location.
#' 
#' During UDDs' first generation, a variety of dose-averaging estimators was developed, with the one proposed by Wetherill et al. (1966) eventually becoming the most popular. This estimator uses only doses observed at *reversal* points: points with a negative response following a positive one, or vice versa. More recent research (Kershaw 1985, 1987; Oron et al. 2022, supplement) strongly indicates that in fact it is better to use all doses starting from some cut-point, rather than skip and choose only reversals. 
#' 
#' The `reversals()` utility identifies reversal points, while `reversmean()` produces a dose-averaging estimate whose starting cut-point is determined by a reversal. User can choose whether to use all doses from that cut-point onwards, or only the reversals as in the older approaches.
#' 
#' More broadly, dose-averaging despite some advantages is not very robust, and also lacks an interval estimate with reliable coverage. Therefore, `reversmean()` provides neither a confidence interval nor a standard erro. Instead, for UDD target estimation we recommend using centered isotonic regression, available via `quickInverse()` in the `cir` package. See Oron et al. 2022 (both article and supplement) for further information.
#' 
#' @references 
#'  - Kershaw CD: A comparison of the estimators of the ED50 in up-and-down experiments. *J Stat Comput Simul* 1987; 27:175–84.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50.
#'  Wetherill GB, Chen H, Vasudeva RB: Sequential estimation of quantal response curves: A new method of estimation. *Biometrika* 1966; 53:439–54

#' 
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x` or shorter by 1, and must be coded `TRUE/FALSE` or 0/1.
#' @param rstart the reversal point from which the averaging begins. See Details.
#' @param all logical: from the starting point onwards, should all values of `x` be used (`TRUE`, default), or only reversal points?
#' @param before logical: whether to start the averaging one step earlier than the starting reversal point. Default `FALSE`, and ignored when `all=FALSE`.
#' @param minfrac a fraction in \eqn{0,1} indicating the minimum fraction of the vector `x` to be used in averaging, in case reversal `rstart` occurs too late in the experiment. Default 0.5.
#' @param full logical: should more detailed information be returned, or only the estimate? (default \code{FALSE})

#' 
#' @return For `reversals()`, the indices of reversal points. For `reversmean()`, if `full=FALSE` returns the point estimate and otherwise returns a data frame with the estimate as well, as the index of the cutoff point used to start the averaging.

#'  
#' @export
#'  
reversmean <- function(x, y, rstart=3, all=TRUE, before=FALSE,
                       maxExclude=0.5,  full=FALSE)
{
# vals
checkDose(x)
checkResponse(y)
n=length(x)
if (!(length(y) %in% c(n-1,n))) stop('X vector must be equal-length or 1 longer than Y.\n')

checkNatural(rstart, toolarge = floor(n/2))
checkTarget(maxExclude, tname = 'maxExclude')
# /vals

revpts=reversals(y)

#### exception handling: 
# if zero reversals, we err out 
if(length(revpts)==0) stop('No reversals. Experiment likely too short, or data-quality error.\n')
# part-degenerate: fewer revs than expected
if(rstart > length(revpts)) rstart=length(revpts) 
# Late start: reverting to some minimal start point:
if(revpts[rstart] > n*maxExclude) revpts[rstart] = floor(n*maxExclude)

# The estimate is anti-climactic:
est=ifelse(all,mean(x[(revpts[rstart] - as.integer(before)):n]),mean(x[revpts[rstart:length(revpts)]]))
if(!full) return(est)
data.frame(est = est, cutoff = revpts[rstart] - as.integer(before) )
}


#' @rdname reversmean
#' @export
#' 
reversals <- function(y) 
{
  which(diff(y)!=0)+1
}




##' Average with adaptive rstarting-point
##' Based on an unpblished concept from Oron (2007)

adaptmean<-function(x,minfrac=2/3,before=FALSE,full=FALSE)
{
# Degenerate case
if(length(unique(x))==1) return(x[1])

n=length(x)
# Means of the tail only; tailmeans[1] is mean of everything, tailmeans[n]=x[n].
tailmeans=rev(cumsum(rev(x))/(1:n))
spacing=mean(diff(sort(unique(x))))
# If you're near where you started, quick exit using the entire sample:
if(abs(tailmeans[2]-x[1])<=spacing) return (tailmeans[1])

signvec=sign(x[-n]-tailmeans[-1])
hinge=min(which(signvec!=signvec[1]))
# Rolling one step backwards, to *before the crossing
if(before) hinge=hinge-1
if(signvec[1]==0) hinge=2 # perfect storm

### Return
if(full) return(list(startpt=hinge,signsmeans=rbind(tailmeans,c(signvec,NA))))
# Applying minimum fraction
minstart=floor(n*(1-minfrac))

ifelse(hinge<minstart,tailmeans[hinge],tailmeans[minstart])
}

### Choi (1971,1990) point+interval estimate 

choiEst<-function(x,y,rstart=2,conf=0.9,full=FALSE)
{
n=length(x)
if (!(length(y) %in% c(n-1,n))) stop('X vector must be equal-length or 1 longer than Y.\n')

revpts=reversals(y)
k=length(revpts)

# The silly 'trough-peak' correction
exes=x[ revpts[rstart:k] ] + 1/2 - y[ revpts[rstart:k] ]
# the ubiquitous k-1:
k1=length(exes)

# Point estimate: these shifted reversals 
pest = mean(exes)

# rho: following Choi (90), we just use the lag-1 ACF to stay out of trouble.
rho0 = acf(exes,lag=1,plot=FALSE)$acf[2]

# Estimating single-obs rho and sigma FWIW (Choi 90 eq. 7-8)

aa = exes[1]^2 + exes[k1]^2
bb = sum(exes[-1]*exes[-k1])
cc = sum(exes[-c(1,k1)]^2)
rhoot = function(x,aa,bb,cc,k)  
	(bb-cc*x)*(1-x^2) - x*(aa + cc*(1+x^2) - 2*bb*x)/k

# eq. 7	
rho = uniroot(rhoot, interval=c(-1,1), aa=aa,bb=bb,cc=cc,k=k1)$root
# eq. 8
sig2 = ( aa + cc*(1+rho^2) - 2*bb*rho ) / k1
sig20 = ( aa + cc*(1+rho0^2) - 2*bb*rho0 ) / k1

# Putting it together (Choi 90 eq. 6)!
eyes = 1:(k1-1)
se = sqrt ( (sig2 / (k1*(1-rho^2) ) ) * (1 + 2*sum( rho^eyes * (k1-eyes) )/k1 ) )
se0 = sqrt ( (sig20 / (k1*(1-rho0^2) ) ) * (1 + 2*sum( rho0^eyes * (k1-eyes) )/k1 ) )

return(c(pest,rho,rho0,k1,se,se0))
}



### Averaging standard error, ad-hoc

avgHalfCI <- function(x,conf=0.9,refq=c(.1,.9),full=FALSE)
{
neff=max(table(x))-1
sdeff=diff(quantile(x,refq,type=6))/2
if(!full) return(qt(0.5+conf/2,df=neff-1)*sdeff/sqrt(neff))
return(data.frame(neff=neff,sdeff=sdeff))
}





