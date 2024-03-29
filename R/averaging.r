### Up-and-Down target estimates based on dose averaging 


#' Original Dixon and Mood (1948) point estimate
#'
#' Basic version; formula assumes uniform spacing but should work anyway
#'
#' In their documentation of the Up-and-Down algorithm, Dixon and Mood (1948) presented an estimation method based on tallying and averaging responses, choosing to use only the positive or negative responses (the less-common of the two), since they reasoned the two mirror each other. This is not strictly true: it ignores both leading/trailing "tail" sequences of identical responses, and repeated visits to the boundary dose in case there are dose boundaries.
#' 
#' The Dixon-Mood estimate (sometimes called Dixon-Massey) is provided here mostly for historical reasons and comparative-simulation uses, and also because this estimate is apparently *(and unfortunately)* still in use in some fields. **It should not be used for actual target-dose estimation in real experiments.** It behaves very poorly even under minor deviations from the most optimal conditions (see, e.g., simulations in the supplement to Oron et al. 2022.).
#' 
#' In order to discourage from actual use in experiments, we do not provide a method for the Dixon-Mood estimator's confidence interval, even though the original article did include one. The interval estimate behaves even more poorly than the point estimate. 
#' 
#' For UDD target estimation we recommend using centered isotonic regression, available via \code{\link{udest}}, an up-and-down adapted wrapper to `cir::quickInverse()`. See Oron et al. 2022 (both article and supplement) for further information, as well as the `cir` package vignette.
#' 
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' 
#' @example inst/examples/avgExamples.r

#' @seealso
#'   - \code{\link{udest}}, the recommended estimation method for up-and-down targets. 
#'   - \code{\link{reversmean}}, The most commonly-used dose-averaging approach (*not* recommended; the recommended one is \code{\link{udest}} referenced above).	 

#' 
#' @export
#' 
#' @return The point estimate

#' @references 
#'  - Dixon WJ, Mood AM. A method for obtaining and analyzing sensitivity data. *J Am Stat Assoc.* 1948;43:109-126.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. [See in particular the open-access Supplement.](https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf)


#' @inheritParams reversmean
#' 
#' @param flip logical: should we flip D-M's approach and use the more-common outcome? (default \code{FALSE})

dixonmood <- function(x, y, full=FALSE, flip=FALSE)
{
n=length(x)
if (!(length(y)==n)) stop('Mismatched lengths.\n')
if(length(unique(y))==1) return(NA) # degenerate case, all 0's or 1's
xvals = sort(unique(x))
spacing = mean(diff(xvals))
# D-M method uses the less frequent outcome, hence might drop leading/trailing "tails"
chosen = ifelse(mean(y)<0.5, 1, 0)
if(flip) chosen = 1-chosen
track = x[y==chosen]

# After all this prep, the formula is anticlimactic :)
if(full) return(c(A=sum((track-min(track))/spacing), N=length(track), d=spacing))

mean(track) + spacing * (0.5-chosen)
}

#----------------------- Utilities for reversal-anchored averaging -----------------#



#' Reversal-anchored averaging estimators for Up-and-Down
#' 
#' Dose-averaging target estimation for Up-and-Down experiments, historically the most popular approach, but not recommended as primary nowadays. Provided for completeness.
#' 
#' Up-and-Down designs (UDDs) allocate doses in a random walk centered nearly symmetrically around a balance point. Therefore, a modified average of allocated doses could be a plausible estimate of the balance point's location.
#' 
#' During UDDs' first generation, a variety of dose-averaging estimators was developed, with the one proposed by Wetherill et al. (1966) eventually becoming the most popular. This estimator uses only doses observed at *reversal* points: points with a negative response following a positive one, or vice versa. More recent research (Kershaw 1985, 1987; Oron et al. 2022, supplement) strongly indicates that in fact it is better to use all doses beginning from some cutoff point, rather than skip most of them and choose only reversals. 
#' 
#' The `reversals()` utility identifies reversal points, whereas `reversmean()` produces a dose-averaging estimate whose cutoff point (which should perhaps be called the *'cut-on'* point) is determined by a reversal. User can choose whether to use all doses from that cut-point onwards, or only the reversals as in the older approaches. A few additional options make the estimate even more flexible.
#' 
#' Starting version 0.2.0, a bootstrap confidence interval (CI) is also provided. See \code{\link{dfboot}} for additional optional parameters to pass to the bootstrap routine, beyond the confidence level `conf`. To skip CI estimation, set `conf = NULL`.   

#' `reversmean()` is compatible mostly with median-targeting UDDs such as the "Classical" (traditional) design of Dixon and Mood. 
#' For general UDD target estimation, particularly off-median targeting designs, we recommend using centered isotonic regression, available via \code{\link{udest}}, an up-and-down adapted wrapper to `cir::quickInverse()`. See Oron et al. 2022 (both article and supplement) for further information, as well as the `cir` package vignette.
#' 
#' @references 
#'  - Kershaw CD: A comparison of the estimators of the ED50 in up-and-down experiments. *J Stat Comput Simul* 1987; 27:175–84.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. [See in particular the open-access Supplement.](https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf)
#'  - Wetherill GB, Chen H, Vasudeva RB: Sequential estimation of quantal response curves: A new method of estimation. *Biometrika* 1966; 53:439–54
#'  
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}
#'
#' @example inst/examples/avgExamples.r

#' @seealso 
#'  - \code{\link{udest}}, the recommended estimation method for up-and-down targets.
#'  - \code{\link{dynamean}}, an unpublished but arguably better approach to dose-averaging (this is *not* the recommended method though; that would be \code{\link{udest}} referenced above).
#' 
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x` or shorter by 1, and must be coded `TRUE/FALSE` or 0/1. `dynamean()` only uses `y` for bootstrap confidence intervals.
#' @param rstart the reversal point from which the averaging begins. Default 3, considered a good compromise between performance and robustness. See Details.
#' @param all logical: from the cutoff point onwards, should all values of `x` be used (`TRUE`, default), or only reversal points as in the Wetherill et al. approach? If set to `FALSE`, then the `before` flag also defaults to `FALSE` regardless of user choice.
#' @param before logical: whether to start the averaging one observation earlier than the cutoff point. Default `FALSE`.
#' @param maxExclude a fraction in \eqn{0,1} indicating the maximum initial fraction of the vector `x` to exclude from averaging, in case the algorithm-identified transition point occurs late in the experiment. Default 1/2.
#' @param full logical: should more detailed information be returned, or only the estimate? (default \code{FALSE})
#' @param conf the CI's confidence level, as a fraction in (0,1). To skip CI calculation set `conf = NULL`.
#' 
#' @return For `reversals()`, the indices of reversal points. For `reversmean()`, if `full=FALSE` returns the point estimate and otherwise returns a data frame with the estimate, as well as the index of the cutoff point used to start the averaging.

#'  
#' @export
#'  
reversmean <- function(x, y, rstart=3, all=TRUE, before=FALSE, conf = 0.9,
                       maxExclude=1/2,  full=FALSE, wetherill = TRUE, ...)
{
# vals
checkDose(x)
checkResponse(y)
n=length(x)
if (!(length(y) %in% c(n-1,n))) stop('X vector must be equal-length or 1 longer than Y.\n')

checkNatural(rstart, toolarge = floor(n/2))
checkTarget(maxExclude, tname = 'maxExclude')
if(!is.null(conf)) checkTarget(conf, tname = 'conf')
# /vals

revpts=reversals(y=y, x=x, directional = wetherill)

#### exception handling: 
# if zero reversals, we return NA (used to just err out) 
if(length(revpts)==0) return(NA)
# stop('No reversals. Experiment likely too short, or data-quality error.\n')

# part-degenerate: fewer revs than expected
if(rstart > length(revpts)) rstart=length(revpts) 
# Late start: reverting to some minimal start point:
if(revpts[rstart] > n*maxExclude) revpts[rstart] = floor(n*maxExclude)

# The estimate is anti-climactic:
est=ifelse(all,mean(x[(revpts[rstart] - as.integer(before)):n]),mean(x[revpts[rstart:length(revpts)]]))

### NEW to 0.2! Bootstrap CI
if(!is.null(conf)) confidence = dfboot(x=x, y=y, conf = conf, full = full, estfun = reversmean,
                              rstart = rstart, all = all, before = before, maxExclude = maxExclude,
                              ...) else confidence = NULL

### Returning
if(full) return(list(est = est, cutoff = revpts[rstart] - as.integer(before), 
                     bootstrap = confidence) )

if(is.null(conf)) return(est)

tmp = c(est, min(confidence[1], est), max(confidence[2], est) ) 
names(tmp) = c('point', names(confidence) )
tmp
}


#' @rdname reversmean
#' @export
#' 
reversals <- function(y, x = NULL, directional = TRUE, ycheck = FALSE) 
{
  if(directional && is.null(x[1])) stop("x must be provided.\n")
  n = length(y)
  if(directional && !(length(x) %in% c(n, n+1) ) ) stop('Mismatched x:y lengths.\n')
     
  if(directional)   
  {
    moves = which(diff(x) != 0)
    hinges = diff(x)[moves]
    tmp = moves[which(diff(hinges)!=0) + 1]
#    cat(tmp,'\n')
#    cat(diff(y[tmp]),'\n')
    if(ycheck) if(any(diff(y[tmp]) == 0)) stop('Inconsistent x:y sequences.\n')
    
  } else tmp = which(diff(y)!=0)+1
  
  tmp
}

#------------------------

#' Up-and-Down target-dose averaging estimate, with dynamic cutoff-point
#'
#' A dose-averaging estimate based on a concept from Oron (2007). Provides a more robust alternative to reversal-based averaging.
#'
#' Historically, most up-and-down studies have used dose-averaging estimates. Many of them focus on reversal points either as anchor/cutoff points -- points where the averaging begins -- or as the **only** doses to use in the estimate.  Excluding doses before the anchor/cutoff is done in order to mitigate the bias due to the arbitrary location of the starting dose. The extent of excluded sample depends on the distance between the starting dose and the 
#'    up-and-down balance point, as well as the random-walk vagaries of an individual experimental run. 

#' 
#' Oron (2007) showed that using only reversals and skipping other doses is generally a bad idea, and also noted that a reversal anchor point is not directly tied to the conceptual motivation for having an anchor/cutoff point.
#'
#'    In practice, some *"lucky"* experiments might not need any exclusion at all (because they started right at the balance point), while others might need to exclude dozens of observations. Reversals do not capture this variability well.
#' 
#' The estimation method coded in `dynamean()` works from a different principle. It identifies **the first crossing point:** the first point at which 
#'   the dose is *"on the other side"* from the starting point, compared with the average of all remaining doses. 
#' The average of all remaining doses is used as a proxy to the (unobservable) balance point. 
#' This approach is far closer to capturing the dynamics described above, and indeed performs well
#'     in comparative simulations (Oron et al. 2022, Supplement).
#'
#' Starting version 0.2.0, a bootstrap confidence interval (CI) is also provided. See \code{\link{dfboot}} for additional optional parameters to pass to the bootstrap routine, beyond the confidence level `conf`. To skip CI estimation, set `conf = NULL`.   

#' The experiment's binary responses (`y`) are only needed as input for confidence interval calculations; otherwise, only the dose-allocation sequence (`x`) is required. 
#' 
#' `dynamean()` is recommended only for median-targeting UDDs, e.g., the "Classical" (traditional) design of Dixon and Mood. For off-median percentiles it might not work as well, and CI coverage will be lacking.
#' 
#' For such targets we recommend using centered isotonic regression, a more robust method available 
#'   together with a confidence interval via \code{\link{udest}}, an up-and-down adapted wrapper to `cir::quickInverse()`.
#'     See Oron et al. 2022 (both article and supplement) for further information, as well as the `cir` package vignette.


#' @inheritParams reversmean

#' @export
#' 
#' @return If `full=FALSE` returns the point estimate. Otherwise returns a list with the index of the cutoff point used to start the averaging, and a 2-row matrix with full sequence of estimates and the tail vs. starting point indicator signs.

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	
#'
#' @example inst/examples/avgExamples.r
#' 
#' @seealso 
#'   - \code{\link{udest}}, the recommended estimation method for up-and-down targets. 
#'   - \code{\link{reversmean}} for the commonly-used reversal-anchored averages mentioned in Details.	  

#' 
#' @references 
#' 
#'  - Oron AP. [*Up-and-Down and the Percentile-finding Problem.*](https://arxiv.org/abs/0808.3004) Ph.D. Dissertation, University of Washington, 2007. 
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. [See in particular the open-access Supplement.](https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf)


dynamean <- function(x, y=NULL, maxExclude=1/2, before=FALSE, full=FALSE, conf = 0.9, ...)
{
  
### Validation
  checkDose(x)
  n=length(x)
if(!is.null(y))
{
  checkResponse(y)
  if (!(length(y) %in% c(n-1,n))) stop('X vector must be equal-length or 1 longer than Y.\n')
}
  checkTarget(maxExclude, tname = 'maxExclude')
  if(!is.null(conf)) checkTarget(conf, tname = 'conf')

  # Degenerate case
if(length(unique(x))==1) return(x[1])

n=length(x)
# Means of the tail only; tailmeans[1] is mean of everything, tailmeans[n]=x[n].
tailmeans=rev(cumsum(rev(x))/(1:n))
spacing=mean(diff(sort(unique(x))))

signvec = sign(x[-n] - tailmeans[-1])
hinge = suppressWarnings( min( which(signvec != signvec[1]) ) )
# Rolling one step backwards, to *before the crossing
if(before) hinge = hinge-1
if(signvec[1]==0) hinge = 2 # perfect conditions

### NEW to 0.2! Bootstrap CI
if(!is.null(conf)) confidence = dfboot(x=x, y=y, conf = conf, 
                    full = full, estfun = dynamean, 
                    before = before, maxExclude = maxExclude, ...) else confidence = NULL

### Return
if(full) return(list(startpt=hinge, signsmeans=rbind(tailmeans,c(signvec,NA)), bootstrap = confidence) )

# Applying minimum fraction
minstart = floor(n * maxExclude)
pest = ifelse(hinge<minstart, tailmeans[hinge], tailmeans[minstart])
if(is.null(conf)) return(pest)

tmp = c(pest, min(confidence[1], pest), max(confidence[2], pest) ) 
names(tmp) = c('point', names(confidence) )
tmp

}





