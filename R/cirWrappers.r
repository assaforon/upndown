#' Centered-Isotonic-Regression (CIR) Estimate for the Up-and-Down Target Dose
#' 
#' Centered Isotonic Regression (CIR) is an extension of isotonic regression (IR), substantially improving upon IR's estimation performance in the dose-response and dose-finding contexts (Oron and Flournoy 2017, Flournoy and Oron 2020). CIR is the recommended method for estimating up-and-down targets. 
#' 
#'  CIR and related methods are available in the `cir` package. The `udest()` function in the present package provides a convenient wrapper for `cir::quickInverse()`, with arguments already set to the appropriate values for estimating the target dose after an up-and-down experiment. The function also returns a confidence interval as default.
#' 
#' **WARNING!** You should not estimate target doses too far removed from the design's actual balance point (definitely no further than 0.1, e.g., estimating the 33rd percentile for a design whose balance point is the median). As Flournoy and Oron (2020) explain, observed response rates are biased away from the balance point. Even though `udest()` performs the rudimentary bias correction described in that article, practically speaking this correction's role is mostly to expand the confidence intervals in response to the bias. It cannot guarantee to provide reliable off-balance-point estimates.
#'
#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' 
#' @export
#'
#' @return A one-row data frame with 4 variables: `target`, `point` (the point estimate), `lowerXYconf, upperXYconf` (the confidence bounds, with `XY` standing for the percents, default `90`).
#'  
#' @param x numeric vector: sequence of administered doses, treatments, stimuli, etc.
#' @param y numeric vector: sequence of observed responses. Must be same length as `x`, and must be coded `TRUE/FALSE` or 0/1.
#' @param target The target response rate for which target dose estimate is requested. Must be a single number in \eqn{(0,1).}
#' @param balancePt In case the design's inherent balance point differs somewhat from `target`, specify it here to improve estimation accuracy. See Details for further explanation. Otherwise, this argument defaults to be equal to `target`.
#' @param conf The desired confidence level for the confidence interval. Default \eqn{90\%.} We do not recommend increasing to \eqn{95\%} unless you have \eqn{\sim 100} or more observations. 
#' @param ... Pass-through argument added for flexible calling context.
#' 
#' @references 
#'  - Oron AP, Flournoy N.  Centered Isotonic Regression: Point and Interval Estimation for Dose-Response Studies. *Statistics in Biopharmaceutical Research* 2017; 9, 258-267. [Author's public version available on arxiv.org.](https://arxiv.org/pdf/1701.05964.pdf)
#'  - Flournoy N, Oron AP. Bias Induced by Adaptive Dose-Finding Designs. *Journal of Applied Statistics* 2020; 47, 2431-2442.  
#'  - Oron AP, Souter MJ, Flournoy N. Understanding Research Methods: Up-and-down Designs for Dose-finding. *Anesthesiology* 2022; 137:137–50. [See in particular the open-access Supplement.](https://cdn-links.lww.com/permalink/aln/c/aln_2022_05_25_oron_aln-d-21-01101_sdc1.pdf)
#'  
#'  @seealso 
#'  \code{\link[cir]{quickInverse}}, `cir` package.
#'  
udest <- function(x, y, target, balancePt = target, conf = 0.9, ...)
{
  require(cir)
  
if(length(target) > 1) stop("Experiment should have a single target.\n")
checkTarget(target)
checkTarget(conf, tname = "'conf'")
if(length(x) != length(y)) stop('x, y, must have same length.\n')
checkDose(x)
checkResponse(y)
if(balancePt[1]!=target)
{
  if(length(balancePt) > 1) stop("Experiment can only have a   single balance point.\n")
  checkTarget(balancePt, tname='balancePt')
  if(abs(balancePt - target) > 0.1) warning("We strongly advise against estimating targets this far from the design's balance point.\n")
}

# And after all this.... it's a one-liner :)

quickInverse(x=x, y=y, target=target, starget = balancePt,
             conf=conf, adaptiveShrink=TRUE,
              adaptiveCurve = (target != 0.5) )

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
#'  `udplot()` is a convenience wrapper to `cir::plot.DRtrace`. This is a base-R plot, so you can use additional options, including preceding the plot command with \code{\link[graphics]{par}} statements. If wishing to save to a file, I recommend utilities such as `png()` or `pdf()`.

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' 
#' @export
#' @import graphics
#' 
#' @seealso 
#'  \code{\link[cir]{plot.DRtrace}}, `cir` package.


#' @inheritParams udest

#' @param shape the plotting shape (DRtrace only): `'circle'` (default), `'square'`, or `'triangle'`.
#' @param connect logical: whether to connect the symbols (generic plotting type `'b'`). Default \code{TRUE} for `udplot()` and \code{FALSE} for `drplot()`.
#' @param symbcol The color of the main plotting symbols and connecting lines. Default 1 (the current palette's first color). Note: if you change the color and inadvertently use \code{col} instead, there might be an error message.
#' @param xtitle,ytitle	x-axis and y-axis titles. Some reasonable defaults are provided, to avoid an annoying error message.
#' @param doselabels (\code{DRtrace} only) Dose values to be plotted along the y-axis. If \code{NULL} (default), those will be the doses in the dataset (i.e., `sort(unique(x))`).
#' @param ...	Other arguments passed on to \code{\link[graphics]{plot}} (e.g., `main` for the main title). 


udplot <- function(x, y, shape='circle', connect=TRUE, symbcol=1, doselabels=NULL, 
                   xtitle = "Observation Order", ytitle = "Dose / Stimulus",...)
{
require(cir)

# val
checkDose(x)
checkResponse(y)

  
plot(DRtrace(x=x, y=y), shape=shape, connect=connect, mcol=symbcol, dosevals=doselabels,
             xlab=xtitle, ylab=ytitle, ...)
  
}

#------------------------------- Dose-Response Plotting -----------------------------

#' Visualizing the dose-response summary of an up-and-down experiment
#'
#'
#'
#'

#' @author Assaf P. Oron \code{<assaf.oron.at.gmail.com>}	  
#' 
#' @export
#' @import graphics
#' 
#' @seealso 
#'  - \code{\link[cir]{plot.doseResponse}}, `cir` package.
#'  - `cir` package vignette.


#' @inheritParams udplot
#' @inheritParams udest

drplot <- function(x, y, shape='X', connect=FALSE, symbcol=1, 
                   addest=FALSE, addcurve=FALSE, target=NULL, balancePt=target, 
                   estcol = 'purple', estsize=2, estsymb=19, esthick=2, curvecol = 'blue',
                   ytitle = "Frequency of Positive Response", xtitle = "Dose / Stimulus",...)
{
require(cir)
  
# val
checkDose(x)
checkResponse(y)
  
  
suppressWarnings(plot(doseResponse(x=x, y=y), pch=shape, connect=connect, mcol=symbcol, 
       xlab=xtitle, ylab=ytitle, ...))

if(addest)
{
  if(is.null(target)) stop("To plot an estimate, please specify the target response rate.\n")
  checkTarget(target)
  est = udest(x=x, y=y, target=target, balancePt=balancePt, ...)
  points(target~point, data=est, pch=estsymb, col=estcol, cex=estsize)
  lines(x = c(est[1,3], est[1,4]), y = rep(target, 2), col=estcol, lwd = esthick)
  
  if(addcurve)
  {
    fest = cirPAVA(x=x, y=y, target=balancePt, adaptiveShrink=TRUE, full=TRUE)
    lines(y ~ x, data = fest$shrinkage, col=curvecol, lwd = esthick)
    
  }
 
}
  
}
