#----- Classical UD Example -----#

# An example used in Oron et al. 2022, Fig. 2.
# It is presented here via the original motivating story:
# "Ketofol"  is a commonly-used anesthesia-inducing mix known to combine its 2 components' 
# beneficial properties, while each component mitigates the other's harmful side-effects. 
# In particular:
#     Propofol reduces blood pressure while ketamine raises it.
# What is *not* known at present, is which mix proportions produce 
# 0 "delta-BP" on average among the population. 

# The classical UD design below administers the mix 0-100% ketamine in 10% increments.
#    The design will concentrate doses around the point where half the population 
#    experiences 0 "delta-BP". (the 'zeroPt' parameter in the code)

doses = seq(0, 100, 10)
m=length(doses) # 11 dose levels

zeroPt=63 # the zero-BP point, in units of percent ketamine

# We assume a Normal ("Probit") dose-response curve,
#   and calculate the value of F (i.e.,  prob (delta BP > 0) at the doses:
equivF = pnorm( (doses - zeroPt) / 20)
round(equivF, 3)

# The vector below represents the values feeding into the Fig. 2B barplot.
# "startdose = 6" means the experiment begins from the 6th out of 11 doses, i.e., a 50:50 mix.

\donttest{
  round(cumulvec(cdf = equivF, matfun = classicmat, startdose = 6, n = 30), 3)
}
# Compare with the *instantaneous* probability distribution to the 30th patient:

round(currentvec(cdf = equivF, matfun = classicmat, startdose = 6, n = 30), 3)
# Classic up-and-down has quasi-periodic behavior with a (quasi-)period of 2. 
# Compare the instantaneous vectors at n=30 and 29:
round(currentvec(cdf = equivF, matfun = classicmat, startdose = 6, n = 29), 3)
# Note the alternating near-zero values. Distributions at even/odd n "communicate"
#    with each other only via the dose boundaries.

# Lastly, the asymptotic/stationary distribution. Notice there is no 'n' argument.

round(pivec(cdf = equivF, matfun = classicmat), 3)

# The cumulative vector at n=30 is not very far from the asymptotic vector. 
# The main difference is that at n=30 there's still a bit more
#    probability weight at the starting dose.
# We can check how much of that extra weight is from the 1st patient, by excluding that data point:

round(cumulvec(cdf = equivF, matfun = classicmat, startdose = 6, n = 30, exclude = 1), 3)


