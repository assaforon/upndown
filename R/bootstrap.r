#' @export
#' 
udboot <- function(x, y, doses =  NULL, estfun = adaptmean, design = krow, desArgs = list(k=1), 
                        target = 0.5, conf = 0.9, B = 1000, seed = NULL, randstart = TRUE,
                        showdots = TRUE, full = FALSE, ...)
{
  require(cir)
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
    bootdat = dfsim(n, starting = startdose, xvals = doses, sprobs = startp, Fvals = bootF, 
                    progress = design, progArgs = desArgs, 
               nlev = length(bootF), ensemble = B, quiet = !showdots, ...)

#### Estimation    
    
    if(identical(estfun, adaptmean)) 
    {
      bootests = apply(bootdat$dose, 2, adaptmean, full=FALSE, conf=NULL, ...)
    } else {
      
      
      
    }
    if(full) return(list(xvals = doses, F = bootF, x = bootdat$dose, ests = bootests) )
    
    tailz = (1-conf)/2
    candout = quantile(bootests, c(tailz, 1-tailz), type = 6)
    names(candout) = paste(c('lower', 'upper'), round(100*conf), 'conf', sep='')
    return(candout)
}