#----- Classical UD Example -----#

# This vector is can be seen as a bar plot in Oron et al. 2022, Fig. 2B.
# It is presented here via the original motivating story:
# "Ketofol" a commonly-used anesthesia-inducing mix
#     known to combine its 2 components' beneficial properties, while each one
#     mitigates the other's harmful side-effects. In particular:
#     Propofol reduces blood pressure while ketamine raises it.
# What is *not* known at present, is which mix proportions produce 0 delta BP
#     on average among the population. 

# The classical UD design below administers the mix 0-100% ketamine in 10% increments
#    The design will concentrate doses around the point where 50% of the population 
#    experiences 0 delta BP. ('zeroPt')

doses = seq(0, 100, 10)
m=length(doses)

zeroPt=63 # percent ketamine
# We assume a Normal ("Probit") dose-response curve,
#   and calculate the value of F (i.e.,  prob (delta BP > 0) at the doses:

equivF = pnorm( (doses - zeroPt) / 20)

# Finally, the values feeding into the Fig. 2B barplot

cumulpi(startdose = 6,equivF,n


