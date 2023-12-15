


avgboot <- function(x, y, doses =  NULL, estfun = adaptmean, target = 0.5, 
                        conf = 0.9, B = 500, dotsout = TRUE, ...)
{
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
    
# F estimate
    cirF = cirPAVA(x = x[1:n], y = y, adaptiveShrink = TRUE, nmin = 1, 
                     target = target, full = TRUE)$output
    bootF = rep(NA, m)
    
    
  
  
  
  
}