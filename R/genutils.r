## GUD target finding
gudtarg<-function(k,u=0,l,b1=1,b2=1)
{
uniroot(f=function(x,kay,you,ell,bee1,bee2) {bee1*pbinom(q=ell,size=kay,prob=x)+bee2*(pbinom(q=you-1,size=kay,prob=x)-1)},interval=0:1,kay=k,you=u,ell=l,bee1=b1,bee2=b2)$root
}

### K-in-a-row targets

k2targ<-function(k,hitarg=TRUE)
{
if(k!=round(k) || k<=0) stop("k can only be a natural number.")
tmp=0.5^(1/k)
if(hitarg) return(tmp)
1-tmp
}

targ2k<-function(targ,expand=1)
{
hi=(targ>=0.5) 
if(!hi) targ=1-targ

kcand=-1/log2(targ)
krange=(floor(kcand-expand):ceiling(kcand+expand))
krange=krange[krange>0]
data.frame(k=krange,balancePt=k2targ(krange,hitarg=hi))
}


############################### Auxiliary utilities


validUDinput<-function(cdf,target)
{
  if(target<=0 || target>=1) stop("Target has to be in (0,1).\n")
  if(min(cdf)<0 || max(cdf)>1 || any(diff(cdf)<0) || var(cdf)==0) stop("cdf should be a CDF.\n")
  
  #ttarg=ifelse(target>0.5,1-target,target)
  if(length(cdf)<3) stop ("These designs don't work with <3 dose levels.\n")
}


