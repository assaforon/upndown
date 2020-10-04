## GUD target finding
gudtarg<-function(k,u=0,l,b1=1,b2=1)
{
uniroot(f=function(x,kay,you,ell,bee1,bee2) {bee1*pbinom(q=ell,size=kay,prob=x)+bee2*(pbinom(q=you-1,size=kay,prob=x)-1)},interval=0:1,kay=k,you=u,ell=l,bee1=b1,bee2=b2)$root
}

