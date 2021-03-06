% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvoutApp.R
\name{conkApp}
\alias{conkApp}
\title{use shiny to explore univariate and multivariate outlier patterns in selected gene sets}
\usage{
conkApp(genesets = conkout::glioSets47)
}
\arguments{
\item{genesets}{a named list of vectors of gene symbols}
}
\description{
use shiny to explore univariate and multivariate outlier patterns in selected gene sets
}
\examples{
if (interactive()) conkApp()
}
