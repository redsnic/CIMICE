% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysis.R
\name{compute_weights_default}
\alias{compute_weights_default}
\title{Compute default weights}
\usage{
compute_weights_default(g, freqs)
}
\arguments{
\item{g}{a graph (must be a DAG with no transitive edges)}

\item{freqs}{observed frequencies of genotypes}
}
\value{
a graph with the computed weights
}
\description{
This procedure computes the weights for edges of a
graph accordingly to CIMICE specification.
(See vignettes for further explainations)
}
\examples{
require(dplyr)
preproc <- example_dataset() \%>\% dataset_preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph_non_transitive_subset_topology(samples, labels)
compute_weights_default(g, freqs)

}
