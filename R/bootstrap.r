#' Generic percentile dose-response / dose-finding bootstrap routine
#' 
#' Bootstrap routine for resampling a dose-finding or dose-response experiment. The bootstrap replicates are generated from a centered-isotonic-regression estimate of the dose-response function, rather than resampled directly. 
#' 
#' 
#' 
#' 
#' 
#' 
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x` or shorter by 1, and must be coded `TRUE/FALSE` or 0/1. `dynamean()` only uses `y` for bootstrap confidence intervals.
#' @param doses the complete set of dose values that *could* have been included in the experiment. Must include all unique values in `x`.

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
#' @seealso \code{\link{simulate}}


#' @export
#' 
udboot <- function(x, y, doses =  NULL, estfun = dynamean, design = krow, desArgs = list(k=1), 
                        target = 0.5, conf = 0.9, B = 1000, seed = NULL, randstart = TRUE,
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
                     target = target)
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
               nlev = length(bootF), ensemble = B, quiet = !showdots) )

    # "Dressing up" the dose levels (which are 1:m in the progress loop above) with real values
    bootdoses = suppressMessages(plyr::mapvalues(bootdat$doses, 1:length(bootF), doses) )
    
    
#### Estimation    
    
    if(identical(estfun, dynamean)) 
    {
      bootests = apply(bootdat$doses, 2, dynamean, full=FALSE, conf=NULL, ...)
    } else {
      
      bootests = rep(NA, B)
      for (a in 1:B) bootests[a] = estfun(x = bootdoses[,a], y = bootdat$responses[,a], 
                        full=FALSE, conf=NULL, target = target, allow1extra = TRUE, ...)
    }
    if(full) return(list(xvals = doses, F = bootF, x = bootdoses, ests = bootests) )
    
    tailz = (1-conf)/2
    candout = quantile(bootests, c(tailz, 1-tailz), type = 6, na.rm = TRUE)
    names(candout) = paste(c('lower', 'upper'), round(100*conf), 'conf', sep='')
    return(candout)
}