### dfsim example

# We provide a "toy example" for how randomized-F simulations might work.
# We have been strong advocated of this simulation approach, especially for performance comparison
#   between designs and between estimators. It is far preferable to the "cherry-picked F" approach.
# Unfortunately, the latter is still more commonly found in dose-finding literature.

# At core is a randomized ensemble of F(x) ("dose-response") curves. It is far easier
#   to generate parametric ensembles rather than nonparametric or semi-parametric ones.
# So this is what we do here. We use the Weibull, being capable of the most diverse ensembles
#   among common parametric families.

# This being a toy example, it is more simplistic than our current practice. For the latter,
#   see either the supplement to Oron et al. 2022, or the NE Journal of SDS due online
#   late 2024 or early 2025.  Ok, let's go!

# We simulate 7-level experiments.
m = 7
# And 10 curves in total
B = 10
# Parameter generation (I chose the parameter bounds after some trial and error)
# Note I am not fixing a seed, so each time you run this you'll get a different ensemble!
wparams = cbind(2^runif(B, -2, 2.5), runif(B, 1, 10) )
round(wparams, 2)
# Now the F(x) curves; they should be in columns
ensemble = apply(wparams, 1, function(x, doses = m) 
                  pweibull(1:doses, shape = x[1], scale = x[2]) )

# Let's see what we got!
round(ensemble, 3)
# The experiment will be "Classical" median-targeting. In real life if 0.5 falls outside
#     of the range of doses, it's not great. For simulation it's fine; it'll enable us to 
#     watch the allocations for such curves heap near the edge with F closest to 0.5.

# Let the experiments be a measly n=15, to keep runtime under the menacing 5 seconds.
# We start all runs smack in the middle, level 4:
sout = dfsim(n = 15, starting = 4, Fvals = ensemble, design = krow, desArgs = list(k=1) )
names(sout)

# "scenarios" are the F values we provided, while "sample" is the set of randomized thresholds
#   simulated from the uniform distribution. If you run comparisons, we recommend that
#   you only randomize the first set in the comparison, and feed the very same thresholds to 
#   all the rest.

# Anyhoo, here are the simulated trajectories. Note that there are n+1 x values in each one, because after 
#     seeing x_n and y_n, the n+1-th allocation can be determined.
sout$doses
# Compare with the F ensemble. Is it what you expected? 
# Probably more random-walky than you thought :)

# The binary responses are below. Before going big on the simulation, it's a good idea to
#    sanity-check and see that the doses and responses match according to the design's rules.
sout$responses
# Some meta details about the design and simulation settings are available here:
sout$details


