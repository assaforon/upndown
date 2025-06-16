## ----setup, include = FALSE------------------------------------------------------------------------------------------------
#### Note to self: Use this before committing to GitHub, after a stable version
# tools::buildVignettes(dir = '.',  tangle=TRUE)
# Then, not sure whether it's necessary to move/copy the files into inst/doc 
# (from https://community.rstudio.com/t/browsevignettes-mypackage-saying-no-vignettes-found/68656/6)

knitr::opts_chunk$set(
  collapse = TRUE, fig.width = 9, fig.height = 7, out.width = 900, out.height = 700, 
  comment = "#"
)

## ----gorla6_1--------------------------------------------------------------------------------------------------------------
library(upndown)

# From Gorla et al. (2017) Table 6
gorla751x = 39 + c(3:0, 1, 2, 1:3, 2, 3, 2, 3)
# With UD data, one can usually discern the y's (except the last one) from the x's:
gorla751y =  c( (1 - diff(gorla751x) ) / 2, 1)

udplot(x=gorla751x, y=gorla751y, main = "Gorla et al. 2017, Material 751", 
       ytitle = "Load (kN)", las = 1)

legend('bottomright', legend = c(expression('Survived 10'^7*' cycles'), 
                              'Failed'), pch=c(1, 19), bty = 'n')


## ----gorla6_2--------------------------------------------------------------------------------------------------------------

drplot(x=gorla751x, y=gorla751y, addest = TRUE, target = 0.5, addcurve = TRUE, 
       percents = TRUE, main = "Gorla et al. 2017, Material 751", 
       xtitle = "Load (kN)", ytitle = "Percent Failure", las = 1)
legend('bottomright', legend = c('CIR target estimate and 90% CI',
              'CIR load-failure curve estimate'), lty = 1, 
              col = c('purple', 'blue'), bty = 'n')
       
udest(gorla751x, gorla751y, target = 0.5)

## ----gorla6_3a-------------------------------------------------------------------------------------------------------------
# Note the manual addition of "dose n+1"
reversmean(c(gorla751x, 41), gorla751y, rstart = 1, all = TRUE, conf = NULL)

## ----gorla6_3b-------------------------------------------------------------------------------------------------------------
# Note the manual addition of "dose n+1"
reversmean(c(gorla751x, 41), gorla751y, rstart = 1, all = TRUE, 
                                   target = 0.5, design = krow, desArgs = list(k=1) )


## ----gorla6_4, width = 100-------------------------------------------------------------------------------------------------
udest(x = gorla751x, y = gorla751y, target = 0.05, balancePt = .5)

## ----gorla9_1, fig.width=13, fig.height=7, out.width=1300, out.height=700, echo = -1---------------------------------------
# Saving current settings as now required by the CRAN powers-that-be :0
op <- par(no.readonly = TRUE)

par(mfrow=1:2, mar=c(4,4,4,1))
# From Gorla et al. (2017) Table 9
gorla951x = 35 + c(1:0, 1:4, 3:2, 3:0, 1, 2, 1)
gorla951y =  c( (1 - diff(gorla951x) ) / 2, 1)

udplot(x=gorla951x, y=gorla951y, main = "Gorla et al. 2017, Material 951", 
       ytitle = "Load (kN)", las = 1)

drplot(x=gorla951x, y=gorla951y, addest = TRUE, target = 0.5, addcurve = TRUE, 
       percents = TRUE, main = "Gorla et al. 2017, Material 951", 
       xtitle = "Load (kN)", ytitle = "Percent Failure", las = 1)

udest(gorla951x, gorla951y, target = 0.5)

par(op) # Back to business as usual ;)

## ----gorla9_2, width = 100-------------------------------------------------------------------------------------------------
dixonmood(x=gorla951x, y=gorla951y)

## ----george1, echo = -1, fig.width = 10, fig.height = 6, out.width=1000, out.height=600------------------------------------
dlabel = 'Phenylephrine dose (micrograms)'
george10x = 80 + 20 * c(1, rep(2, 5), 1, 1, 0, 0, rep(1, 7), 0:2, 2, 2, rep(1, 4), 2, 1, 1, 2, 2, 
                        rep(3, 5), 4, 5, 5, rep(4, 6))
george10y = c(ifelse(diff(george10x) > 0, 0, 1), 1)
udplot(x=george10x, y=george10y, ytitle = dlabel )


## ----george2, echo = -1, fig.width = 12, fig.height = 7, out.width=1300, out.height=700------------------------------------
drplot(x=george10x, y=george10y, addest = TRUE, target = 0.9, addcurve = TRUE, balancePt = 10/11, xtitle = dlabel)

