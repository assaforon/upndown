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


### Let us plot the dose-allocation time series.

# Saving current settings as now required by the CRAN powers-that-be :0
op <- par(no.readonly = TRUE)

par(mar=c(4,4,4,1), mgp=c(2.5,0.8,0), cex.axis = 0.7, las = 1)
udplot(doses, responses, main='Van Elstraete et al. 2008 Study', 
       xtitle = "Patient Number", ytitle = 'Gabapentin (mg/kg)') 


#' Overlay the ED50 reported in the article (21.7 mg/kg):
abline(h = 21.7)

#' The authors cite a little-known 1991 article by Dixon as the method source.
#' However, in their author rejoinder they claim to have used the Dixon-Mood (1948) estimate.


# Our package does include the Dixon-Mood point estimate.
#  (w/o the CIs, because we do not endorse this estimation approach)
# Does it reproduce the article estimate?
dixonmood(doses, responses)

# Not at all! Let us overlay this one in red
abline(h = dixonmood(doses, responses), col=2)

# We have found that many articles claiming to use Dixon-Mood (or Dixon-Massey) actually
# Do something else. For example, in this article they report that 
#   "it is necessary to reject sequences with three to six identical results".
# Nothing like this appears in the original Dixon-Mood article, where the estimation method
#   involves identifying the less-common response (either 0 or 1), and using only x values
#   associated with these responses; obviating the need to exclude specific sequences.
#
# More generally, these historical estimates have long passed their expiry dates. 
#   Their foundation is not nearly as solid as, e.g., linear regression, 
#      and it's time to stop using them.

# That said, our package does offer two more types of dose-averaging estimates.
# Both are able to take advantage of the "n+1" dose-allocation, which is determined by
#    the last dose and response:
n = length(doses)
dosePlus1 = doses[n] + ifelse(responses[n]==0, 1, -1)
reversmean(c(doses, dosePlus1), responses)
# Interestingly, in this particular case the answer is very similar to the Dixon-Mood estimate.

# The `reversmean()` default averages all doses from the 3rd reversal point onwards.
# By the way, at what point did the third reversal happen? 
#     It'll be the 3rd number in this vector:
reversals(responses)

# Far more commonly in literature, particularly in sensory studies, 
#   one encounters the 1960s-era approach (led by Wetherill) of taking *only doses  
#   at reversal points, usually starting from the first one. `reversmean()` can do that too:
wetherill = reversmean(c(doses, dosePlus1), responses, all = FALSE, rstart = 1)
wetherill
# This one gives an even lower result than the previous ones.
abline(h = wetherill, col = 3)

# There's another approach to dose-averaging, although it is not in use anywhere that we know of.
# It does not require the y values at all. The underlying assumption is that the dose 
#   sequence has done enough meandering around the true balance point, to provide information
#   about where (approximately) the starting-dose effect is neutralized.
adaptmean(c(doses, dosePlus1))
# Again a bit curiously, this relatively recent approach gives a result similar to what
#   the authors reported (but not similar to the original Dixon-Mood).
# This is not too surprising, since here `adaptmean()` excludes the first one-third of doses,
#   which is approximately what happened if indeed the authors excluded all those long dose-increase
#   sequences at the start.

# All this shows how dicey dose-averaging, at face value a simple and effective method, can become.
# The sample size here is rather large for up-and-down studies, and yet because of the unlucky
#    choice of starting point (which in many studies, due to safety concerns cannot be evaded)
#    there is really no good option of which observations to exclude.

# This is one reason why we strongly recommend using Centered Isotonic Regression as default:
defest = udest(doses, responses, target = 0.5)
abline(h = defest$point, col = 'purple')
# For this dataset, it is the highest of all the estimates.

legend('bottomright', col = c(1:3, 'purple'), 
       legend = c("Article's estimate", 'Dixon-Mood', 'Reversals (Wetherill)', 'Standard (CIR)'), 
       lty = 1, bty='n', cex = 0.8)

par(op) # Back to business as usual ;)
