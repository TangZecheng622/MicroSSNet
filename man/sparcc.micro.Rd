% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreatMatrix.R
\name{sparcc.micro}
\alias{sparcc.micro}
\title{SparCC Microbial Correlation}
\usage{
sparcc.micro(data, R = 10, ncpus = 1)
}
\arguments{
\item{data}{Numeric matrix, input data with features in columns and samples in rows.}

\item{R}{Integer, number of bootstrap replicates.}

\item{ncpus}{Integer, number of CPU cores to use.}
}
\value{
A list containing the correlation matrix and the p-value matrix.
}
\description{
This function calculates the SparCC correlation matrix and p-values.
}
