## Identify reversals in an UD experiment's time series
reversals<-function(y) which(diff(y)!=0)+1

##' Original Dixon-Mood 48 point estimate
##' Basic version; formula assumes uniform spacing but should work anyway
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

#' Reversal-anchored averaging estimators: reversal-only and reversal-startpoint
reversmean<-function(x,y,start=3,all=TRUE,before=FALSE,full=FALSE)
{
n=length(x)
if (!(length(y) %in% c(n-1,n))) stop('X vector must be equal-length or 1 longer than Y.\n')

revpts=reversals(y)
# exception handling
if(length(revpts)==0) return(mean(x[-1]))
if(start>length(revpts)) start=length(revpts)
# The est is anti-climactic:
est=ifelse(all,mean(x[(revpts[start]-before):n]),mean(x[revpts[start:length(revpts)]]))
if(!full) return(est)
data.frame(est=est,cutoff=revpts[start]-before)
}

##' Average with adaptive starting-point
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
# Return
if(full) return(list(startpt=hinge,signsmeans=rbind(tailmeans,c(signvec,NA))))
# Applying minimum fraction
minstart=floor(n*(1-minfrac))
ifelse(hinge<minstart,tailmeans[hinge],tailmeans[minstart])
}
