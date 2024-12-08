#' Generic percentile dose-response / dose-finding bootstrap routine
#' 
#' Bootstrap routine for resampling a dose-finding or dose-response experiment. The bootstrap replicates are generated from a centered-isotonic-regression (CIR) estimate of the dose-response function, rather than resampled directly. 
#' 
#' The function should be able to generates bootstrap resamples of any dose-finding design, as long as `design, desArgs` are specified correctly.  See \code{\link{dfsim}} for the requirements from such a function.
#' 
#' Like Chao and Fuh (2001) and Stylianou et al. (2003), the bootstrap samples are generated indirectly, by estimating a dose-response curve F from the data, then generating an ensemble of bootstrap experiments using the same design used in the original experiment. Unlike these two which used parametric or isotonic regression, respectively, with no bias-mitigation and no additional provisions to improve coverage, our implementation uses centered isotonic regression with bias-mitigation 
#' 
#' 
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x` or shorter by 1, and must be coded `TRUE/FALSE` or 0/1. 
#' @param doses the complete set of dose values that *could* have been included in the experiment. Must include all unique values in `x`.
#' @param estfun the estimation function to be bootstrapped. Default \code{\link{dynamean}}
#' @param design,desArgs design details passed on to \code{\link{dfsim}}; the former is a function and the latter a list of its arguments and values. For self-consistent bootstrapping, this must specify the design used in the actual experiment. See \code{\link{dfsim}}. 
#' @param target The target percentile to be estimated (as a fraction). Again must be the same one estimated in the actual experiment. Default 0.5.
#' @param balancePt In case the design's inherent balance point differs somewhat from `target`, specify it here to improve estimation accuracy. See Details for further explanation. Otherwise, this argument defaults to be equal to `target`.
#' @param B Size of bootstrap ensemble, default 1000.
#' @param seed Random seed; default `NULL` which leads to a "floating" seed, varying between calls.
#' @param randstart Logical: should the bootstrap runs randomize the starting dose, or use the same starting dose as the actual experiment? Default `TRUE`, which we expect to produce better properties. The randomization will be weighted by the real data's dose-specific sample sizes.
#' @param showdots Logical: should "progress dots" be printed out as the bootstrap runs progress? Default `TRUE`
#' @param full Logical: controls how detailed the output is. Default (`FALSE`) is only the resulting interval bounds, while `TRUE` returns a list with the full bootstrap ensemble of doses, responses and estimates, as well as the generating dose-response curve and the bootstrap's dose set.
#' @param conf the CI's confidence level, as a fraction in (0,1).
#' @param ... Additional parameters passed on to estimation functions.
#' 
#' 
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
#' @seealso \code{\link{dfsim}}


#' @export
#' 
dfboot <- function(x, y, doses =  NULL, estfun = dynamean, design, desArgs, 
                        target, balancePt = target, conf = 0.9, B = 1000, seed = NULL, randstart = TRUE,
                        showdots = TRUE, full = FALSE, ...)
{
  requireNamespace('cir')
  requireNamespace('plyr')
  
    ### validation
  
    if( !is.null(doses) & !all(x %in% doses) ) stop("'doses' must include all x values in the experiment.\n")
    checkResponse(y)
    n=length(y)
    if (!(length(x) %in% c(n, n+1))) stop('X vector must be equal-length or 1 longer than Y.\n')
    
### Setting up dose set for shrinkage of F  

    dose0 = sort(unique(x))
    
    if(is.null(doses)) # No user input; creating a "double sandwiched" dose set
    {
      spaces = diff(dose0)
      doses = c( min(dose0)-2*spaces[1], min(dose0)-spaces[1], 
                dose0, max(dose0)+rev(spaces)[1], max(dose0)+2*rev(spaces)[1] )
      
    } else doses = sort(unique(doses)) 
# Getting full-set indices of used doses  
    m = length(doses)
    indices = match(dose0, doses)
    if(any(is.na(indices))) stop('"doses" does not include all dose values in x!\n')
    minused = min(indices)
    maxused = max(indices)
    
### F estimate for the bootstrapping
    cirF = cir::cirPAVA(x = x[1:n], y = y, adaptiveShrink = TRUE, nmin = 1, 
                     target = balancePt)
    bootF = rep(NA, m)
    bootF[indices] = cirF

    if(minused > 1) bootF[(minused-1)] = bootF[minused] / 2
    if(minused > 2) bootF[1:(minused-2)] = 0
    if(maxused < m) bootF[(maxused+1)] = (1 + bootF[maxused]) / 2
    if(maxused < m-1) bootF[(maxused+2):m] = 1
    
    if(any(is.na(bootF))) bootF[is.na(bootF)] = 
        approx(doses[!is.na(bootF)], bootF[!is.na(bootF)], xout = doses[is.na(bootF)] )$y
    
#    return(bootF)

#### Calling dfsim() to generate ensemble
    startdose = NULL
    if(!randstart) {
      startdose =  match(x[1], doses)
    } else {
      tmp = table(x)
      startp = rep(0, length(doses))
      startp[match(names(tmp), doses)] = tmp/sum(tmp)
    }
    bootdat = suppressMessages( dfsim(n, starting = startdose, sprobs = startp, Fvals = bootF, 
                    design = design, desArgs = desArgs, 
                    ensemble = B, showdots = showdots) )

    # "Dressing up" the dose levels (which are 1:m in the progress loop above) with real values
    bootdoses = suppressMessages(plyr::mapvalues(bootdat$doses, 1:length(bootF), doses) )
    
    
#### Estimation    
    
    if(identical(estfun, dynamean)) 
    {
      bootests = apply(bootdoses, 2, dynamean, full=FALSE, conf=NULL, ...)
    } else {
      
      bootests = rep(NA, B)
      for (a in 1:B) bootests[a] = estfun(x = bootdoses[,a], y = bootdat$responses[,a], 
                        full=FALSE, conf=NULL, target = target, allow1extra = TRUE, ...)
    }
    if(full) return(list(xvals = doses, F = bootF, x = bootdoses, y = bootdat$responses, ests = bootests) )
    
    tailz = (1-conf)/2
    candout = quantile(bootests, c(tailz, 1-tailz), type = 6, na.rm = TRUE)
    names(candout) = paste(c('lower', 'upper'), round(100*conf), 'conf', sep='')
    return(candout)
}