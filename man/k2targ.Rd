% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genutils.r
\name{k2targ}
\alias{k2targ}
\alias{ktargOptions}
\alias{g2targ}
\title{Up-and-Down Target Calculation and Calibration}
\usage{
k2targ(k, hitarg = TRUE)

ktargOptions(target, tolerance = 0.1)

g2targ(cohort, lower, upper)
}
\description{
Up-and-down target calculation and design options/guidance given a user-desired target.
}
\details{
This suite of utilities helps users
\itemize{
\item Figure out the approximate target response-rate given design parameters
\item Suggest or specify design parameters, given user's target response-rate.
}

Up-and-down designs (UDDs) generate random walks over dose space, with most dose-allocations usually taking place near the design's de-facto target percentile, called the \strong{"balance point"} by some theorists to distinguish it from the user's designated target (Oron and Hoff 2009, Oron et al. 2022).

Most k-in-a-row and group UDD parameter combinations yield balance points that are irrational percentiles of the dose-response function, and therefore are unappealing as official experimental targets.

However, since the UD dose distribution has some width, and since even the balance point itself is only a close approximation for the actual average of allocated doses, the user's target \strong{does not have to be identical to the balance point.} It only needs to be \emph{"close enough"}.

The \code{k2targ()} and \code{g2targ()} utilities are intended for users who already have a specific k-in-a-row or group design in mind, and only want to verify its balance point. The complementary utilities \verb{ktargOptions(), gtargOptions()} provide a broader survey of design-parameter options within user-specified constraints, given a desired target.

Lastly, \code{bcoin()} returns the biased-coin probabilities given the user's designated target. In contrast to the two other UDDs described above, the biased-coin design can target any percentile with a precisely matched balance point. That said, k-in-a-row and group UDDs offer some advantages over biased-coin in terms of properties and operational simplicity.
}