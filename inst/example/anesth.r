### Some UD experiment in anesthesiology

## Examples 1a-b: from
## Columb MO, Lyons, G. 1995. Determination of the Minimum Local Analgesic Concentrations of Epidural Bupivacaine and Lidocaine in Labor, Anesthesia & Analgesia 81, 833-837.
## open access: https://journals.lww.com/anesthesia-analgesia/Fulltext/1995/10000/Determination_of_the_Minimum_Local_Analgesic.30.aspx

## Creating dose-response objects using the classes from package 'cir'

# The experimental trace for bupivacaine 
# doses and responses are in the order given
# (doses are in multiples of 0.01%w/v)

cl95bupi=DRtrace(x=0.01*c(15:9,10:6,7,6,7:5,6,7:4,5,6,5:8,7,6),
	y=c(rep(c(rep(1,5),0),2),1,0,1,1,0,0,rep(1,3),0,0,1,rep(0,3),1,1,0))

layout(matrix(1:4,nrow=2,byrow=TRUE),widths=3:2)
par(mar=c(4,4,4,1),mgp=c(2.5,.8,0))
plot(cl95bupi,main='Columb and Lyons 95: Bupivacaine') 
# Compare with the article's Figure 1; the commented line below makes it look more similar
# plot(cl95bupi,shape='square',connect=FALSE,cex=2)

cl95bupiRates=doseResponse(cl95bupi)
plot(cl95bupiRates,main='Bupivacaine Dose-Response')
# Due to the adaptive design, these observed rates are actually biased. 
# See the rates after an empirical bias correction:
points(DRshrink(cl95bupiRates,target=0.5),col=2)

# The experimental trace for lidocaine 
# doses and responses are in the order given
# (doses are in multiples of 0.05%w/v)
# Note the first few doses use double spacing

cl95lido=DRtrace(x=0.05*c(seq(20,12,-2),11:9,10:7,8,7,8,7:9,8:6,7,8,7,8:5,6,7),
	y=c(rep(1,7),0,rep(1,3),0,1,0,1,0,0,rep(1,3),0,0,1,0,rep(1,3),rep(0,3)))
	
plot(cl95lido,main='Columb and Lyons 95: lidocaine') 
# Compare with the article's Figure 2; the commented line below makes it look more similar
# plot(cl95lido,shape='square',connect=FALSE,cex=2)

cl95lidoRates=doseResponse(cl95lido)
plot(cl95lidoRates,main='Lidocaine Dose-Response')
# Due to the adaptive design, these observed rates are actually biased. 
# See the rates after an empirical bias correction:
points(DRshrink(cl95lidoRates,target=.5),col=2)

##################### Articles from Pace and Stylianou

p03ropi=DRtrace(x=0.01*c(10:8,9,10:7,8,9:7,8,7:11,10:7,8,7,8,7:12,11,12:10),
	y=c(1,1,0,0,1,1,1,0,0,1,1,0,1,rep(0:1,each=4),0,1,0,1,rep(0,5),1,0,1,1,1))
p03ropiRates=doseResponse(p03ropi)






