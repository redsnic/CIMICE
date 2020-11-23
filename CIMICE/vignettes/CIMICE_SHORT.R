## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(CIMICE)
library(dplyr)
library(igraph)

## ----cars---------------------------------------------------------------------
# read from file
head(read.CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE)))

## -----------------------------------------------------------------------------
# from a string
read.CAPRI.string("
s\\g A B C D
S1 0 0 0 1
S2 1 0 0 0
S3 1 0 0 0
S4 1 0 0 1
S5 1 1 0 1
S6 1 1 0 1
S7 1 0 1 1
S8 1 1 0 1
")

## -----------------------------------------------------------------------------
# using CIMICE::make.dataset and CIMICE::update_df
# genes
make.dataset(A,B,C,D) %>% 
    # samples
    update_df("S1", 0, 0, 0, 1) %>%
    update_df("S2", 1, 0, 0, 0) %>%
    update_df("S3", 1, 0, 0, 0) %>%
    update_df("S4", 1, 0, 0, 1) %>%
    update_df("S5", 1, 1, 0, 1) %>%
    update_df("S6", 1, 1, 0, 1) %>%
    update_df("S7", 1, 0, 1, 1) %>%
    update_df("S8", 1, 1, 0, 1)

## -----------------------------------------------------------------------------
preproc <- dataset.preprocessing(example.dataset())
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]

## -----------------------------------------------------------------------------
g <- graph.non.transitive.subset.topology(samples,labels)

## ---- echo=FALSE--------------------------------------------------------------
V(g)$vertex.size <- rep(10, length(V(g)))
plot(g, vertex.size=rep(55, length(V(g))))

## -----------------------------------------------------------------------------
W <- compute.weights.default(g, freqs)

## -----------------------------------------------------------------------------
draw.visNetwork(g, W, labels)

## -----------------------------------------------------------------------------
quick.run(example.dataset())

