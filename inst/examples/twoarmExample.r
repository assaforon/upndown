
# Data from: Camorcia et al. (2011), Eur. J Anaesth., 28, 240–244.

# ED50 of anesthetic for routine surgery: men vs. women
c11m_x = c(4:6, 5, 4:7, 6, 5, 6:4, 5:8, 7, 8 , 7, 8, 7:9, 8:11, 10, 9)
c11m_y = c( (1 - sign(diff(c11m_x)) )/ 2, 0) 

c11f_x = c(4, 5, 4:6, 5:3, 4:7, rep(6:5, 3), 6, 7:4, 5, 6:4, 5, 6, 5)
c11f_y = c( (1 - sign(diff(c11f_x)) )/ 2, 1) 

#### The confidence interval (equivalent to a significance test)
diff11ci95_mf = armDiffCI(c11m_x, c11m_y, c11f_x, c11f_y, target = 0.5)	
diff11ci95_mf

# NOTE: the experiment actually had 3 arms, so officially some multiple-testing adjustment is required.
# However, as shown above the contrast we tested was non-significant at 95%, even without such an adjustment.
# Unlike our re-analysis, the original authors using a very optimistic CI method have this contrast 
#    to be significant, even at the alpha = 0.05/3 level.

#### Optional: Data visualization

layout(matrix(1:4, nrow = 2), widths = 3:2)
par(mar = c(4,4,4,2), mgp = c(2.6, 0.5, 0), tck = -0.01, las = 1, 
    cex.lab = 1.8, cex.axis = 1.4, cex.main = 2)

dosename = "Bupivacaine (mg)"
dosevals = 2:11
arms = c("Women",  "Men")
types = c('Trajectory', 'Dose-Response')

# Overlap of the two 84% CIs is roughly equivalent to a 95% test for the CI of the difference.
# So we use this in visualizing each arm's uncertainty
ovconf = 0.84


udplot(x = c11f_x, y = c11f_y, main = paste(arms[1],types[1]), ytitle = dosename, 
       doselabels = dosevals, ylim = range(dosevals), cex = 2 )
udplot(x = c11m_x, y = c11m_y, main = paste(arms[2],types[1]), ytitle = dosename, 
       doselabels = dosevals, ylim = range(dosevals), cex = 2 )
legend('topleft',legend=c('Effective','Ineffective'),pch=c(19,1),bty='n', cex = 1.4)


drplot(x = c11f_x, y = c11f_y, main = paste(arms[1],types[2]), xtitle = dosename, 
       ytitle = 'Proportion Effective',  obsize = 1., curvethick = 2,
       doselabels = dosevals, addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5)
drplot(x = c11m_x, y = c11m_y, main = paste(arms[2],types[2]), xtitle = dosename, 
       ytitle = 'Proportion Effective', obsize = 1., curvethick = 2,
       doselabels = dosevals, addest = TRUE, addcurve = TRUE, conf = ovconf, target = 0.5)

