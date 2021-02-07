## GUD target finding
gudtarg<-function(k,u=0,l,b1=1,b2=1)
{
uniroot(f=function(x,kay,you,ell,bee1,bee2) {bee1*pbinom(q=ell,size=kay,prob=x)+bee2*(pbinom(q=you-1,size=kay,prob=x)-1)},interval=0:1,kay=k,you=u,ell=l,bee1=b1,bee2=b2)$root
}

### K-in-a-row targets

k2targ<-function(k,hiend=TRUE)
{
if(k!=round(k) || k<=0) stop("k can only be a natural number.")
tmp=0.5^(1/k)
if(hiend) return(tmp)
1-tmp
}

targ2k<-function(targ,expand=1)
{
hi=(targ>=0.5) 
if(!hi) targ=1-targ

kcand=-1/log2(targ)
krange=(floor(kcand-expand):ceiling(kcand+expand))
krange=krange[krange>0]
data.frame(k=krange,balancePt=k2targ(krange,hiend=hi))
}

