
#'  **An up-and-down experiment that has generated some controversy**
#'  
#' Van Elstraete, AC et al. The Median Effective Dose of Preemptive Gabapentin 
#'      on Postoperative Morphine Consumption After Posterior Lumbar Spinal Fusion. 
#'      *Anesthesia & Analgesia* 2008, 106: 305-308.

# It was a classical median-finding up-and-down study.

doses = c(4:7, 6:13, 12:19, 18:21, 20, 19:23, 22, 21:23, 22:19, 20:23, 
          22:24, 23, 22, 23, 22:25, 24:22, rep(23:24,2), 23, 22)
# With U&D, responses (except the last one) can be read off the doses:
responses = c( (1 - sign(diff(doses)))/2, 0 )


#' ### Plots plots plots!

layout(t(1:2), widths=3:2)
par(mar=c(4,4,4,1), mgp=c(2.5,0.8,0), cex.axis = 0.7, las = 1)

#' The experimental trajectory / time-series / "trace" (pick your favorite name!)
#' Note the changed argument names for x and y axis titles
udplot(doses, responses, main='', 
        xtitle = "Patient Number", ytitle = 'Gabapentin (mg/kg)') 
#' Compare with the article's Figure 1; the commented line below makes it look more similar
# udplot(doses, responses, shape='square', connect=TRUE)

# The dose-response plot, rarely encountered in U&D articles. 
# We can also add the CIR estimate right there:
drplot(doses, responses, main=' Dose-Response', percents = TRUE,
       addest = TRUE, target = 0.5, addcurve = TRUE,
       xtitle = 'Gabapentin (mg/kg)', ytitle = "Percent Effective")

#' ### Estimates

#' Let us actually see the numbers of those Centered-Isotonic-Regression (CIR) estimates!
#' Note that our default confidence-interval is 90%. Change it via the 'conf' argument.

udest(doses, responses, target = 0.5)
#' Compare with the article: 21.7 mg/kg (95% CI 19.9–23.5). 
#' They cite a little-known 1991 article by Dixon as the method source.
#' However, in their author rejoinder they claim to have used the Dixon-Mood estimate.
#' 
#' ## Toy example of plotting a group UD dataset
#' 
#' Also showing off some udplot() options
#' 
#' Not an actual experiment (made-up data)
#' The design is purportedly GUD (3,0,1), targeting the 20th percentile
#' 

gsize = 3
x = rep(c(1:3, 2:4), each = gsize)
y = c(rep(0, 8), 1, rep(0,7), 1, 1)

udplot(x=x, y=y, cohort=gsize, connect=FALSE, shape='triangle')
