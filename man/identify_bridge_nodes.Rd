% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BridgeNet.R
\name{identify_bridge_nodes}
\alias{identify_bridge_nodes}
\title{Identify Bridge Nodes}
\usage{
identify_bridge_nodes(cor_matrix, table1name, table2name)
}
\arguments{
\item{cor_matrix}{Matrix, correlation matrix between the two domains.}

\item{table1name}{Character, name prefix for species in the first domain.}

\item{table2name}{Character, name prefix for species in the second domain.}
}
\value{
A vector of node names that act as bridges between the two domains.
}
\description{
Extracts nodes directly related between two domains based on the correlation matrix.
}
