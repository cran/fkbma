% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s3-methods.R
\name{plot.rjMCMC}
\alias{plot.rjMCMC}
\title{Plotting function for \code{rjMCMC} results}
\usage{
\method{plot}{rjMCMC}(
  x,
  ...,
  variables = NULL,
  sample_type = "fitted",
  effect_type = "exposure_effect",
  plot_type = "cred",
  level = 0.95,
  aux_vars = list(),
  facet_by = NULL,
  pip_cutoff = 0.1
)
}
\arguments{
\item{x}{A fitted model object from \code{rjMCMC}.}

\item{...}{Additional arguments to be passed to other methods or functions.}

\item{variables}{A vector of variable names to include in the plot. Automatically set to continuous variables if \code{NULL}.}

\item{sample_type}{Character string specifying the type of sample: "fitted", "predictive", or "estimand".
"fitted" and "predictive" are compatible with \code{plot_type = "cred"}.
"estimand" is compatible with \code{plot_type = "hist"} or "trace" (only used for individual parameter trajectories).}

\item{effect_type}{Character string indicating the effect type: "exposure_effect" or "outcome".
For "exposure_effect", the function plots the fitted or predictive effect of exposure; for "outcome", it plots
outcome values without interaction with exposure.}

\item{plot_type}{Character string specifying the plot type: "cred" for credible interval plots, or "hist"/"trace" for histogram or trace plots
of individual parameters (only for \code{sample_type = "estimand"}).}

\item{level}{Numeric value for the credible interval level (default is 0.95).}

\item{aux_vars}{A list of auxiliary variables and their fixed values. Each element name must match a model variable.}

\item{facet_by}{A vector of variable names to facet by in the plot. Automatically set to binary model variables if \code{NULL}.}

\item{pip_cutoff}{Numeric threshold for the posterior inclusion probability (PIP) of model variables to include in the plot.}
}
\value{
A \code{ggplot2} object or a grid of plots.
}
\description{
This function generates plots for model results from \code{rjMCMC} based on specified sample type, effect type, and plot type.
The function is flexible for various combinations of \code{sample_type}, \code{effect_type}, and \code{plot_type}, as outlined below.
}
\details{
\itemize{
\item \strong{Sample and Plot Compatibility}:
\itemize{
\item For \code{sample_type = "estimand"}, only \code{plot_type = "hist"} or "trace" is allowed, as these are designed to visualize the posterior distribution
or MCMC trajectory of individual parameters. Parameters like \code{intercept}, \code{trt}, and \code{sigma} are agnostic to \code{effect_type} as they
do not interact with exposure.
\item \code{plot_type = "cred"} is designed for use with \code{sample_type = "fitted"} or "predictive" and shows credible intervals for the outcome (y-axis)
across biomarker values (x-axis) by covariate pattern. \code{effect_type} controls whether the exposure effect or main effect is displayed.
}
\item \strong{Effect Types}:
\itemize{
\item \code{outcome} plots either the fitted or predictive values without exposure interaction, allowing for exposure (\code{trt}) values to be specified.
\item \code{exposure_effect} plots the interaction of the exposure effect across different covariate patterns.
}
}
}
\examples{
\donttest{
# Example dataset
data("simulated_data")

candsplinevars <- c("X_1")
candbinaryvars <- paste0("Z_", 1:5)
candinter <- c(candsplinevars, candbinaryvars)

results <- rjMCMC(simulated_data, candsplinevars, candbinaryvars, candinter,
                  outcome = "Y", factor_var = "trt")
plot(results, sample_type = "fitted", effect_type = "exposure_effect", plot_type = "cred")
plot(results, sample_type = "estimand", plot_type = "hist")
}
}
