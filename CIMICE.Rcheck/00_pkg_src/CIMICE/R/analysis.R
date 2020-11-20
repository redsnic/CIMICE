# This script contains functions to quicken the
# data analysis and to prepare working examples
utils::globalVariables(c("A","B","C","D","freq"))

#' Creates a simple example dataset
#'
#' @return a simple mutational matrix
#'
#' @examples
#' example.dataset()
#'
#' @export example.dataset
example.dataset <- function(){
    make.dataset(A,B,C,D) %>% # genes
        # samples
        update_df("S1", 0, 0, 0, 1) %>%
        update_df("S2", 1, 0, 0, 0) %>%
        update_df("S3", 1, 0, 0, 0) %>%
        update_df("S4", 1, 0, 0, 1) %>%
        update_df("S5", 1, 1, 0, 1) %>%
        update_df("S6", 1, 1, 0, 1) %>%
        update_df("S7", 1, 0, 1, 1) %>%
        update_df("S8", 1, 1, 0, 1)
}

#' Creates a simple example dataset with frequency column
#'
#' @return a simple mutational matrix
#'
#' @examples
#' example.dataset.withFreqs()
#'
#' @export example.dataset.withFreqs
example.dataset.withFreqs <- function(){
    dataset <- make.dataset(A,B,C,D) %>% # genes
        # samples
        update_df("G1", 0, 0, 0, 1) %>%
        update_df("G2", 1, 0, 0, 0) %>%
        update_df("G3", 1, 0, 0, 1) %>%
        update_df("G4", 1, 1, 0, 1) %>%
        update_df("G5", 1, 0, 1, 1) %>%
        update_df("G6", 1, 1, 0, 1)
    freq <- c(1,2,1,2,1,1)
    dataset %>% cbind(freq)
}

#' Run CIMICE preprocessing
#'
#' executes the preprocessing steps of CIMICE
#'
#' Preprocessing steps:
#'
#' 1) dataset is compacted
#'
#' 2) genotype frequencies are computed
#'
#' 3) labels are prepared
#'
#' @param dataset a mutational matrix as a dataframe
#'
#' @return a list containing the mutational matrix ("samples"),
#' the mutational frequencies of the genotypes ("freqs"),
#' the node labels ("labels") and finally the gene names ("genes")
#'
#' @examples
#' require(dplyr)
#' example.dataset() %>% dataset.preprocessing
#'
#' @export dataset.preprocessing
dataset.preprocessing <- function(dataset){
    # compact
    compactedDataset <- compact.dataset.easy(dataset)
    # remove frequencies column from the compacted dataset
    samples <- as.matrix(compactedDataset %>% select(-freq))
    # save genes' names
    genes = colnames(samples)
    # keep the information on frequencies for further analysis
    freqs = as.matrix(compactedDataset %>% ungroup() %>% select(freq))
    freqs = freqs/sum(freqs)
    freqs = c(freqs,0)
    # prepare node labels listing the mutated genes for each node
    labels <- prepare.labels(samples, genes)
    # fix Colonal genotype absence, if needed
    fix <- fix.clonal.genotype(samples, freqs, labels)
    samples = fix[["samples"]]
    freqs = fix[["freqs"]]
    labels = fix[["labels"]]
    # return a list with the prepared dataset and its additional information
    list("samples" = samples, "freqs" = freqs,
        "labels" = labels, "genes" = genes)
}

#' Run CIMICE preprocessing for poulation format dataset
#'
#' executes the preprocessing steps of CIMICE
#'
#' Preprocessing steps:
#'
#' 1) genotype frequencies are computed
#'
#' 2) labels are prepared
#'
#' @param dataset a mutational matrix as a dataframe (with freq column)
#'
#' @return a list containing the mutational matrix ("samples"),
#' the mutational frequencies of the genotypes ("freqs"),
#' the node labels ("labels") and finally the gene names ("genes")
#'
#' @examples
#' require(dplyr)
#' example.dataset.withFreqs() %>% dataset.preprocessing.population
#'
#' @export dataset.preprocessing.population
dataset.preprocessing.population <- function(dataset){
    # dataset is already compacted per hypothesys
    compactedDataset <- dataset
    # remove frequencies column from the compacted dataset
    samples <- as.matrix(compactedDataset %>% select(-freq))
    # save genes' names
    genes = colnames(samples)
    # keep the information on frequencies for further analysis
    freqs = as.matrix(compactedDataset %>% ungroup() %>% select(freq))
    freqs = freqs/sum(freqs)
    freqs = c(freqs,0)
    # prepare node labels listing the mutated genes for each node
    labels <- prepare.labels(samples, genes)
    # fix Colonal genotype absence, if needed
    fix <- fix.clonal.genotype(samples, freqs, labels)
    samples = fix[["samples"]]
    freqs = fix[["freqs"]]
    labels = fix[["labels"]]
    # return a list with the prepared dataset and its additional information
    list("samples" = samples, "freqs" = freqs,
        "labels" = labels, "genes" = genes)
}

#' Default preparation of graph topology
#'
#' By default, CIMICE computes the relation between
#' genotypes using the subset relation.
#' For the following steps it is also important
#' that the transitive edges are removed.
#'
#' @param samples mutational matrix
#' @param labels genotype labels
#'
#' @return a graph with the wanted topology
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' graph.non.transitive.subset.topology(samples, labels)
#' @export graph.non.transitive.subset.topology
graph.non.transitive.subset.topology <- function(samples, labels){
    # compute edges based on subset relation
    edges <- build.topology.subset(samples)
    # remove transitive edges and prepare igraph object
    build.subset.graph(edges, labels)
}

#' Compute default weights
#'
#' This procedure computes the weights for edges of a
#' graph accordingly to CIMICE specification.
#' (See vignettes for further explainations)
#'
#' @param g a graph (must be a DAG with no transitive edges)
#' @param freqs observed frequencies of genotypes
#'
#' @return a graph with the computed weights
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' compute.weights.default(g, freqs)
#'
#' @export compute.weights.default
compute.weights.default <- function(g, freqs){
    # prepare adj matrix
    A <- as.matrix(as_adj(g))
    # pre-compute exiting edges from each node
    no.of.children <- get.no.of.children(A,g)
    # compute the four steps
    upWeights <- computeUPW(g, freqs, no.of.children, A)
    normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
    downWeights <- computeDWNW(g, freqs, no.of.children, A, normUpWeights)
    normDownWeights <- normalizeDWNW(g, freqs, no.of.children, A, downWeights)
    normDownWeights
}

#' Run CIMICE defaults
#'
#' This function executes CIMICE analysis on a dataset using default settings.
#'
#' @param dataset a mutational matrix as a data frame
#'
#' @return a visNetwork object representing the graph computed by CIMICE
#'
#' @examples
#' quick.run(example.dataset())
#'
#' @export quick.run
quick.run <- function(dataset){
    # preprocess
    preproc <- dataset.preprocessing(example.dataset())
    samples <- preproc[["samples"]]
    freqs   <- preproc[["freqs"]]
    labels  <- preproc[["labels"]]
    genes   <- preproc[["genes"]]
    g <- graph.non.transitive.subset.topology(samples,labels)
    W <- compute.weights.default(g, freqs)
    draw.visNetwork(g, W, labels)
}
