% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{computeDWNW}
\alias{computeDWNW}
\title{Down weights computation}
\usage{
computeDWNW(g, freqs, no.of.children, A, normUpWeights)
}
\arguments{
\item{g}{graph (a Directed Acyclic Graph)}

\item{freqs}{observed genotype frequencies}

\item{no.of.children}{number of children for each node}

\item{A}{adjacency matrix of G}

\item{normUpWeights}{normalized up weights as computed by normalizeUPW}
}
\value{
a vector containing the Up weights for each edge
}
\description{
Computes the Down weights formula using
a Dinamic Programming approach (starting call),
see vignettes for further explaination.
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
# prepare adj matrix
A <- as.matrix(as_adj(g))
# pre-compute exiting edges from each node
no.of.children <- get_no_of_children(A,g)
upWeights <- computeUPW(g, freqs, no.of.children, A)
normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
computeDWNW(g, freqs, no.of.children, A, normUpWeights)

}
