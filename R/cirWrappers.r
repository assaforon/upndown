#' Centered-Isotonic-Regression (CIR) Estimate for the Up-and-Down Target Dose
#' 
#' Centered Isotonic Regression (CIR) is an extension of isotonic regression (IR), substantially improving upon IR's estimation performance in the dose-response and dose-finding contexts (Oron and Flournoy 2017, Flournoy and Oron 2019). CIR is the recommended method for estimating up-and-down targets. 
#' 
#'  CIR and related methods are available in the `cir` package. The `udest()` function in the present package provides a convenient wrapper for `cir::quickInverse()`, with arguments already set to the appropriate values for estimating the target dose after an up-and-down experiment. The function also returns a confidence interval as default.
#' 
#' **WARNING!** You should not estimate target doses too far removed from the design's actual balance point (definitely no further than 0.1, e.g., estimating the 33rd percentile for a design whose balance point is the median). As Flournoy and Oron (2019) explain, observed response rates are biased away from the balance point. Even though `udest()` performs the rudimentary bias correction described in that article, practically speaking this correction's role is mostly to expand the confidence intervals in response to the bias. It cannot guarantee to provide reliable off-balance-point estimates.
#'
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	 
#' 
#' @example inst/examples/wrapExamples.r

#' 
#' @export
#'
#' @return Geneally, a one-row data frame with 4 variables: `target`, `point` (the point estimate), `lowerXYconf, upperXYconf` (the confidence bounds, with `XY` standing for the percents, default `90`).
#' 
#' However, if `conf = NULL` only the point estimate will be returned. This is for compatibility with bootstrap confidence intervals (which we use for dose-averaging estimators but not for CIR), and with UD ensemble simulation in general.
#'  
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x`, and must be coded `TRUE/FALSE` or 0/1.
#' @param target The target response rate for which target dose estimate is requested. Must be a single number in \eqn{(0,1).}
#' @param balancePt In case the design's inherent balance point differs somewhat from `target`, specify it here to improve estimation accuracy. See Details for further explanation. Otherwise, this argument defaults to be equal to `target`.
#' @param conf The desired confidence level for the confidence interval (CI). Default \eqn{90\%.} We do not recommend increasing to \eqn{95\%} unless you have \eqn{\sim 100} or more observations. Setting to `NULL` if you don't want a CI. See more under "Value".
#' @param curvedCI logical: should confidence-interval boundaries rely upon an outwardly-curving interpolation (`TRUE`) or linear? If `NULL` (default), it will be `TRUE` for targets outside the 40th-60th percentile range.
#' @param allow1extra logical: allow `length(x)` to be either equal or 1 greater than `length(y)`? (default `FALSE`) The *"n+1"* dose-allocation, determined from the last allocations and responses, might be tagged onto `x`. If this point is provided and `allow1extra=TRUE`, `udplot()` will show it as a grey diamond; the other functions will ignore it.
#' @param ... Pass-through argument added for flexible calling context.
#' 
#' @references 
#'  - Oron AP, Flournoy N.  Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. *Statistics in Biopharmaceutical Research* 2017; 9, 258-267. [Author's public version available on arxiv.org.](https://arxiv.org/pdf/1701.05964)
#'  - Flournoy N, Oron AP. Bias Induced by Adaptive Dose-Finding Designs. *Journal of Applied Statistics* 2019; 47, 2431-2442.  
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. [See in particular the open-access Supplement.](https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf)

#'  
#' @seealso 
#'  - \code{\link[cir]{quickInverse}}, `cir` package.
#'  - `cir` package vignette.

#### The function:

udest <- function(x, y, target, balancePt = target, conf = 0.9, allow1extra = FALSE, curvedCI = NULL, ...)
{
  requireNamespace('cir')

# validation
  
if(length(target) > 1) stop("Experiment should have a single target.\n")
checkTarget(target)
if(!is.null(conf)) checkTarget(conf, tname = "'conf'")
n = length(y)
if(length(x) != n) 
{
  if (length(x) == n+1 && allow1extra) x = x[1:n] else
  stop('x, y, must have same length.\n')
}
checkDose(x)
checkResponse(y)
if(balancePt[1]!=target)
{
  if(length(balancePt) > 1) stop("Experiment can only have a single balance point.\n")
  checkTarget(balancePt, tname='balancePt')
  if(abs(balancePt - target) > 0.1) warning("\nWe strongly advise against estimating targets this far from the design's balance point.\n")
}

if(is.null(conf)) confi = 0.9 else confi = conf

# And after all this.... it's a one-liner :)
if(is.null(curvedCI)) curvedCI = (target > 0.6 || target < 0.4)

tmp = cir::quickInverse(x=x, y=y, target=target, starget = balancePt,
          conf = confi, adaptiveShrink = TRUE, adaptiveCurve = curvedCI, ... )

if(is.null(conf)) return(tmp$point)

tmp

}

#------------------------------- Trace Plotting -----------------------------

#' Visualizing the time series of an up-and-down experiment
#' 
#' Plotting function for the "trace" (time series) of an up-and-down experiment, showing the observation order on the x-axis, and the dose *(treatment, stimulus, etc.)* strength on the y-axis. Uses utilities from the `cir` package.
#' 
#' This simple and handy visualization approach was presented already by Dixon and Mood (1948). 
#'  - It conveys directly the meaning of *"up-and-down"*, because the administered dose/stimulus strength is on the y-axis, whereas observation order is on the x-axis. 
#'  - Filled symbols stand for positive responses and open symbols for negative.   
#'  - The design's transition rules can be usually inferred directly from the plot.
#'  
#'  `udplot()` is a convenience wrapper to `cir::plot.DRtrace`. This is a base-R plot, so you can use additional options, including preceding the plot command with \code{\link[graphics]{par}} statements, or following up with \code{\link[graphics]{legend}}. When wishing to save to a file, I recommend utilities such as `png()` or `pdf()`.

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' 
#' @example inst/examples/wrapExamples.r

#' @return Returns invisibly after plotting. If you would like to save the plot to a file, embed the plotting code in a standard R graphics export code sequence, (e.g., `pdf(...)` before the plotting function, and `dev.off()` after it).
#' 
#' @export
#' @import graphics
#' 
#' @seealso 
#'  - \code{\link[cir]{plot.DRtrace}}, `cir` package.
#'  - \code{\link{drplot}} for the up-and-down dose-response and estimate plotting.
#'  - `cir` package vignette.
#'  
#' @references 
#'  - Dixon WJ, Mood AM. A method for obtaining and analyzing sensitivity data. *J Am Stat Assoc.* 1948;43:109-126.
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. 



#' @inheritParams udest

#' @param cohort for a group/cohort UD design, the cohort/group size (a single number). In case of variable cohort size, this can be a vector the same length as `x, y`, with each observation's cohort assignment.
#' @param shape the plotting shape (DRtrace only): `'circle'` (default), `'square'`, or `'triangle'`.
#' @param connect logical: whether to connect the symbols (generic plotting type `'b'`). Default \code{TRUE} for `udplot()` and \code{FALSE} for `drplot()`.
#' @param symbcol The color of the main plotting symbols and connecting lines. Default 1 (the current palette's first color). Note: if you change the color and inadvertently use \code{col} instead, there might be an error message.
#' @param xtitle,ytitle	x-axis and y-axis titles. Some reasonable defaults are provided, to avoid an annoying error message.
#' @param doselabels Dose values to be plotted along the y-axis (`udplot()`) or x-axis (`drplot()`). Must be on same scale as `x`. If \code{NULL} (default), those will be the doses in the dataset (i.e., `sort(unique(x))`).
#' @param ...	Other arguments passed on to \code{\link[graphics]{plot}} (e.g., `main` for the main title). 


udplot <- function(x, y, cohort=NULL, shape='circle', connect=TRUE, symbcol=1, 
                   doselabels=NULL, allow1extra = FALSE, 
                   xtitle = "Observation Order", ytitle = "Dose / Stimulus",...)
{

# val
checkDose(x)
checkResponse(y)
n = length(y)
if(length(x) != n) 
{
  if (length(x) != n+1 || !allow1extra) 
    stop('x, y, must have same length.\n')
}

if(!is.null(cohort)) checkNatural(cohort, parname = 'cohort')
xbounds = c(1, ifelse(allow1extra, length(x), n) )
plot(cir::DRtrace(x=x[1:n], y=y, cohort=cohort), shape=shape, connect=connect, mcol=symbcol, dosevals=doselabels,
             xlab=xtitle, ylab=ytitle, xlim = xbounds, ...)
if(length(x)==n+1) points(n+1, x[n+1], pch = 23, bg = 'grey', col = 'grey')
  
}

#------------------------------- Dose-Response Plotting -----------------------------

#' Visualizing the dose-response summary of an up-and-down experiment
#'
#' Dose-response plotting function for up-and-down data, with doses/stimuli on the x-axis, and the proportion of positive responses on the y-axis. Includes an option to plot the target-dose estimate. Uses utilities from the `cir` package.
#' 
#' After an up-and-down experiment, it is highly recommended to plot not just the experiment's *"trace"*  time-series (\code{\link{udplot}}), but also the dose-response summaries. This utility provides a convenient interface for doing that.
#'  - It summarizes the response rates at each participating dose, and plots them. At default, symbol area is proportional to the number of observations at each dose.
#'  - Optionally, the centered-isotonic-regression (CIR) target-dose estimate and its confidence interval are also calculated and plotted.   
#'  - A further option allows for plotting the entire estimated CIR dose-response curve.
#'  
#'  `drplot()` is a convenience wrapper to `cir::plot.doseResponse`, with the added option of plotting the estimate. Some specific options, such as disabling the proportional-area symbol plotting, are accessible only via `plot.doseResponse` arguments (specified in your `drplot()` call and passed through the `...`). 
#'  
#'  This is a base-R plot, so you can use additional options, including preceding the plot command with \code{\link[graphics]{par}} statements, or following up with \code{\link[graphics]{legend}}. When wishing to save to a file, I recommend utilities such as `png()` or `pdf()`.

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' 
#' @export
#' @import graphics
#' 
#' @return Returns invisibly after plotting. If you would like to save the plot to a file, embed the plotting code in a standard R graphics export code sequence, (e.g., `pdf(...)` before the plotting function, and `dev.off()` after it).

#' 
#' @seealso 
#'  - \code{\link[cir]{plot.doseResponse}}, `cir` package.
#'  - \code{\link{udplot}} for the "trace" time-series plot.
#'  - `cir` package vignette.

#' @references 
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. 


#' @inheritParams udplot
#' @inheritParams udest
#' 
#' @param addest logical: should we add the CIR target-dose estimate and its confidence interval? If `FALSE` (default), then arguments `addcurve, target, balancePt, conf, estcol, estsize, estsymb, esthick, curvecol` - are all ignored.
#' @param addcurve logical: should we add the complete estimated CIR dose-response curve? Default `FALSE`, and only relevant when `addest = TRUE`.
#' @param estcol,estsize,estsymb,esthick,curvethick,curvetype,curvecol graphical parameters controlling the colors, symbol choice, size, thickness, of the target-dose and CIR-curve visuals.
#' @param percents logical, whether to represent the y-axis as percents rather than a fraction.

drplot <- function(x, y, shape='X', connect=FALSE, symbcol=1, percents = FALSE, 
                   doselabels=NULL, obsize = NULL,
                   addest=FALSE, addcurve=FALSE, target=NULL, balancePt=target, conf=0.9,
                   estcol = 'purple', estsize=2.5, estsymb=19, esthick=2.5, 
                   curvethick=1.5, curvetype=1, curvecol = 'blue', allow1extra = FALSE,
                   ytitle = "Frequency of Positive Response", xtitle = "Dose / Stimulus",...)
{
requireNamespace('cir')
  
# val
  if(length(target) > 1) stop("Experiment should have a single target.\n")
  checkTarget(target)
  if(!is.null(conf)) checkTarget(conf, tname = "'conf'")
  n = length(y)
  if(length(x) != n) 
  {
    if (length(x) == n+1 && allow1extra) x = x[1:n] else
      stop('x, y, must have same length.\n')
  }
  checkDose(x)
  checkResponse(y)
  
tmp = cir::doseResponse(x=x, y=y)
if(percents) tmp$y = 100 * tmp$y
  
plot(tmp, pch=shape, connect=connect, mcol=symbcol, refsize = obsize,
       xlab=xtitle, ylab=ytitle, dosevals = doselabels, ...)

if(addest)
{
  if(is.null(target)) stop("To plot an estimate, please specify the target response rate.\n")
  checkTarget(target)
  est = udest(x=x, y=y, target=target, balancePt=balancePt, conf=conf, ...)
  if(percents) est$target = 100 * est$target
  points(target ~ point, data=est, pch=estsymb, col=estcol, cex=estsize)
  lines(x = c(est[1,3], est[1,4]), y = rep(est$target, 2), col=estcol, lwd = esthick)
  
  if(addcurve)
  {
    fest = cir::cirPAVA(x=x, y=y, target=balancePt, adaptiveShrink=TRUE, full=TRUE)
    if(percents) fest$shrinkage$y = 100 * fest$shrinkage$y
    lines(y ~ x, data = fest$shrinkage, col=curvecol, lwd = curvethick, lty = curvetype)
    
  }
 
}
  
}
