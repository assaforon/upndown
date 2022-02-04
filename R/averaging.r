### Up-and-Down target estimates based on dose averaging 


##' Original Dixon and Mood (1948) point estimate
#'
##' Basic version; formula assumes uniform spacing but should work anyway
#'
#' In their documentation of the Up-and-Down algorithm, Dixon and Mood (1948) presented an estimation method based on tallying responses


#' @param x time series of doses given
#' @param y time series of binary responses, should be coded 0/1 or \code{FALSE/TRUE}
#' @param full logical: should more detailed information be returned? (default \code{FALSE})
#' @param flip logical: should we flip D-M's approach and use the more-common outcome? (default \code{FALSE})

dixonmood<-function(x,y,full=FALSE,flip=FALSE)
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

##' Identify reversals in an UD experiment's time series
#' @param y time series of binary responses, should be coded 0/1 or \code{FALSE/TRUE}

reversals <- function(y) which(diff(y)!=0)+1

xtrace <- function(x) 
{
trans=diff(x)
c(1,which(trans!=0) + 1)
}


#' Reversal-anchored averaging estimators: reversal-only and reversal-startpoint
reversmean<-function(x,y,rstart=3,all=TRUE,before=0,full=FALSE)
{
n=length(x)
if (!(length(y) %in% c(n-1,n))) stop('X vector must be equal-length or 1 longer than Y.\n')

revpts=reversals(y)
#### exception handling
if(length(revpts)==0) { # fully degenerate, no reversals
	if(full) return(data.frame(est=mean(x[-1]),cutoff=1))
	return(mean(x[-1]))
}
# part-degenerate: fewer revs than expected
if(rstart>length(revpts)) rstart=length(revpts) 
# Exception handling for 'before' (so we don't start before observation 1)
if (revpts[rstart]<=before) before=revpts[rstart]-1

# The estimate is anti-climactic:
est=ifelse(all,mean(x[(revpts[rstart]-before):n]),mean(x[revpts[rstart:length(revpts)]]))
if(!full) return(est)
data.frame(est=est,cutoff=revpts[rstart]-before)
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

# Estimating single-obs sigma FWIW (Choi 90 eq. 8)

aa = exes[1]^2 + exes[k1]^2
bb = sum(exes[-1]*exes[-k1])
cc = sum(exes[-c(1,k1)]^2)
rhoot = function(x,aa,bb,cc,k)  
	(bb-cc*x)*(1-x^2) - x*(aa + cc*(1+x^2) - 2*bb*x)/k
	
rho = uniroot(rhoot, interval=c(-1,1), aa=aa,bb=bb,cc=cc,k=k1)$root
sig2 = ( aa + cc*(1+rho^2) - 2*bb*rho ) / k1

# Putting it together (Choi 90 eq. 6)!
eyes = 1:(k1-1)
se = sqrt ( (sig2 / (k1*(1-rho^2) ) ) * (1 + 2*sum( rho^eyes * (k1-eyes) )/k1 ) )

return(c(pest,rho,rho0,k1,se))
}



### Averaging standard error, ad-hoc

avgHalfCI <- function(x,conf=0.9,refq=c(.1,.9),full=FALSE)
{
neff=max(table(x))-1
sdeff=diff(quantile(x,refq,type=6))/2
if(!full) return(qt(0.5+conf/2,df=neff-1)*sdeff/sqrt(neff))
return(data.frame(neff=neff,sdeff=sdeff))
}





