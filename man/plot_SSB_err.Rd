% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diagnostics.R
\name{plot_SSB_err}
\alias{plot_SSB_err}
\title{Relative Error Plot}
\usage{
plot_SSB_err(output, type, ...)
}
\arguments{
\item{output}{The output matrix from a model run.}

\item{type}{Switch. Type of SSB comparison. Options: "initial", "current", "status".}
}
\value{
Returns a boxplot of relative error.
}
\description{
Plot relative error statistic between SSBs
}
\details{
This function uses data from the \code{output} matrix to plot the relative error between operating and assessment model output quantities.
This will need adaptation as more results are acquired. The main one will be to have the option to vectorize the boxplots like Pauls.
}
