
## This example is from one of the earliest applications of up-and-down in anesthesiology:
## Columb MO, Lyons, G. 1995. Determination of the Minimum Local Analgesic Concentrations of Epidural Bupivacaine and Lidocaine in Labor, Anesthesia & Analgesia 81, 833-837.
## open access: https://journals.lww.com/anesthesia-analgesia/Fulltext/1995/10000/Determination_of_the_Minimum_Local_Analgesic.30.aspx

# It was a classical median-finding up-and-down study.

# Bupivacaine: 
bupix = 0.01 * c(15:9,10:6,7,6,7:5,6,7:4,5,6,5:8,7,6)
# With classical U&D, responses (except the last one) can be read off the doses:
bupiy = c( (1 - sign(diff(bupix)))/2, 0 )

# Lidocaine (note that spacing was halved after observation 5):
lidox = 0.05 * c(seq(20,12,-2),11:9,10:7,8,7,8,7:9,8:6,7,8,7,8:5,6,7)
lidoy = c( (1 - sign(diff(lidox)))/2, 0 )

### Plots plots plots!

layout(matrix(1:4, nrow=2, byrow=TRUE), widths=3:2)
par(mar=c(4,4,4,1), mgp=c(2.5,0.8,0), cex.axis = 0.7, las = 1)

# The experimental trajectory / time-series / "trace" (pick your favorite name!)
# Note the changed argument names for x and y axis titles
udplot(bupix, bupiy, main='Columb and Lyons 95: Bupivacaine', 
        xtitle = "Patient Number", ytitle = 'Bupivacaine (% w/v)') 
# Compare with the article's Figure 1; the commented line below makes it look more similar
# udplot(bupix, bupiy, shape='square', connect=FALSE, cex=2)

# The dose-response plot, rarely encountered in U&D articles. 
#   (Columb and Lyons do provide the info as a table, sort-of - see Table 2)
# We can also add the CIR estimate right there:
drplot(bupix, bupiy, main='Bupivacaine Dose-Response', percents = TRUE,
       addest = TRUE, target = 0.5, addcurve = TRUE,
       xtitle = 'Bupivacaine (% w/v)', ytitle = "Percent Effective")

# Same two plots for lidocaine

udplot(lidox, lidoy, main='Columb and Lyons 95: Lidocaine',
  xtitle = "Patient Number", ytitle = 'Lidocaine (% w/v)') 

drplot(lidox, lidoy, main='Lidocaine Dose-Response', percents = TRUE,
       addest = TRUE, target = 0.5, addcurve = TRUE,
       xtitle = 'Bupivacaine (% w/v)', ytitle = "Percent Effective")

#----------- Let us actually see the numbers of those Centered-Isotonic-Regression (CIR) estimates!
# Note that our default confidence-interval is 90%. Change it via the 'conf' argument.

udest(bupix, bupiy, target = 0.5)
# Compare with the article: 0.065% (95% CI 0.045-0.085) via "Dixon-Massey" (presumably, Dixon-Mood), 
#                       and 0.062% (95% CI 0.052-0.072) via logistic regression

udest(lidox, lidoy, target = 0.5)
# Compare with the article: 0.37% (95% CI 0.32-0.42) via "Dixon-Massey"  
#                       and 0.36% (95% CI 0.31-0.41) via logistic regression

# Our package does include the Dixon-Mood point estimate.
#  (w/o the CIs, because we do not endorse this estimation approach)
# Does it reproduce the article numbers?
dixonmood(bupix, bupiy)
dixonmood(lidox, lidoy)
# Not quite. Perhaps indeed the Dixon-Massey textbook (published a few years after Dixon and Mood)
#    included some tweak on the original. Or perhaps, the formula is a bit open to variations in
#    translation when coding it (originally it was a pencil-and-paper affair).


