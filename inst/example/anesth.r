### Some UD experiment in anesthesiology

## Examples 1a-b: from
## Columb MO, Lyons, G. 1995. Determination of the Minimum Local Analgesic Concentrations of Epidural Bupivacaine and Lidocaine in Labor, Anesthesia & Analgesia 81, 833-837.
## open access: https://journals.lww.com/anesthesia-analgesia/Fulltext/1995/10000/Determination_of_the_Minimum_Local_Analgesic.30.aspx

## Creating dose-response objects using the classes from package 'cir'

# The experimental trace - doses and responses in the order given
# (doses are in multiples of 0.01%w/v)

cl95bupi=DRtrace(x=c(15:9,10:6,7,6,7:5,6,7:4,5,6,5:8,7,6),
	y=c(rep(c(rep(1,5),0),2),1,0,1,1,0,0,rep(1,3),0,0,1,rep(0,3),1,1,0))
layout(t(1:2))
plot(cl95bupi,main='Columb and Lyons 95: Bupivacaine') 
# Compare with the article's Figure 1; the commented line below makes it look more similar
# plot(cl95bupi,shape='square',connect=FALSE,cex=2)

cl95bupiRates=doseResponse(cl95bupi)
plot(cl95bupiRates,main='Bupivacaine Dose-Response')
