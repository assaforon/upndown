#  Let's use an 8-dose design, and  a somewhat asymmetric CDF

exampleF = pweibull(1:8, shape = 2, scale = 4)
# You can plot if you want: plot(exampleF)

# Here's how the transition matrix looks for the median-finding classic up-and-down

round(classicmat(exampleF), 2)
# Note how the only nonzero diagonals are at the opposite corners. That's how 
#   odd-n and even-n distributions communicate (see examples for vector functions).
# Also note how "up" probabilities (the 1st upper off-diagnoal) are decreasing, 
#   while "down" probabilities (1st lower off-diagonal) are increasing, and 
#   start exceeding "up" moves at row 4.

# Now, let's use the same F to target the 90th percentile, which is often
#    the goal of anesthesiology dose-finding studies.
#    We use the biased-coin design (BCD) presented by Durham and Flournoy (1994):

round(bcdmat(exampleF, target = 0.9), 2)

# Note that now there's plenty of probability mass on the diagonal (i.e., repeating same dose).

# Another option, actually with somewhat better operational characteristics, 
#   is "k-in-a-row". Let's see what k to use:

ktargOptions(.9, tolerance = 0.05)

# Even though nominally k=7's target is closest to 0.9, it's generally preferable
#    to choose a somewhat smaller k. So let's go with k=6.
# We must also specify whether this is a low (<0.5) or high target.

round(kmatMarg(exampleF, k = 6, lowTarget = FALSE), 2)

# Compare and contrast with the BCD matrix above! At what dose do the "up" and "down"
#   probabilities flip? 

# Lastly, if you want to see a 43 x 43 matrix - the full state matrix for k-in-a-row, 
#      run the following line:


\donttest{
  round(kmatFull(exampleF, k = 6, lowTarget = FALSE), 2)
}



