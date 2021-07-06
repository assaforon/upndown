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


## Identify reversals in an UD experiment's time series
reversals<-function(y) which(diff(y)!=0)+1


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
