% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weights.R
\name{normalizeDWNW}
\alias{normalizeDWNW}
\title{Down weights normalization}
\usage{
normalizeDWNW(g, freqs, no.of.children, A, downWeights)
}
\arguments{
\item{g}{graph (a Directed Acyclic Graph)}

\item{freqs}{observed genotype frequencies}

\item{no.of.children}{number of children for each node}

\item{A}{adjacency matrix of G}

\item{downWeights}{Down weights as computed by computeDWNW}
}
\value{
a vector containing the normalized Down weights for each edge
}
\description{
Normalizes Down weights so that the sum
of weights of edges exiting a node is 1
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
downWeights <- computeDWNW(g, freqs, no.of.children, A, normUpWeights)
normalizeUPW(g, freqs, no.of.children, A, downWeights)

}
