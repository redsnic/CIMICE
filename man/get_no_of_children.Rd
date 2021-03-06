% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{get_no_of_children}
\alias{get_no_of_children}
\title{Get number of children}
\usage{
get_no_of_children(A, g)
}
\arguments{
\item{A}{Adjacency matrix of the graph g}

\item{g}{a graph}
}
\value{
a vector containing the number of children for each node in g
}
\description{
Compute number of children for each node given an adj matrix
}
\examples{
require(dplyr)
require(igraph)
preproc <- example_dataset() \%>\% dataset_preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph_non_transitive_subset_topology(samples, labels)
A <- as_adj(g)
get_no_of_children(A, g)

}
