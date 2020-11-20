## ---- echo=F------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align="center")

## ---- error=F, message=F, results = "hide"------------------------------------
library(CIMICE)

## ---- error=F, message=F, results = "hide"------------------------------------
# Dataframe manipulation
library(dplyr) 
# Plot display
library(ggplot2)
# Improved string operations
library(glue)
# Dataframe manipulation
library(tidyr)
# Graph data management
library(igraph)
# Remove transitive edges on o graph
library(relations)
# Interactive graph visualization
library(networkD3)
# Interactive graph visualization
library(visNetwork)
# Correlation plot visualization
library(ggcorrplot)
# Functional R programming
library(purrr)
# Graph Visualization
library(ggraph)

## -----------------------------------------------------------------------------
# Read input dataset in CAPRI/CAPRESE format
dataset.big <- read.CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))

## ---- echo=F------------------------------------------------------------------
head(dataset.big)

## -----------------------------------------------------------------------------
# genes
dataset <- make.dataset(A,B,C,D) %>%
    # samples
    update_df("S1", 0, 0, 0, 1) %>%
    update_df("S2", 1, 0, 0, 0) %>%
    update_df("S3", 1, 0, 0, 0) %>%
    update_df("S4", 1, 0, 0, 1) %>%
    update_df("S5", 1, 1, 0, 1) %>%
    update_df("S6", 1, 1, 0, 1) %>%
    update_df("S7", 1, 0, 1, 1) %>%
    update_df("S8", 1, 1, 0, 1) 

## ---- echo=F------------------------------------------------------------------
dataset
as.matrix(dataset)

## -----------------------------------------------------------------------------
gene.mutations.hist(dataset.big)

## -----------------------------------------------------------------------------
sample.mutations.hist(dataset.big, binwidth = 10)

## -----------------------------------------------------------------------------
select.genes.on.mutations(dataset.big, 100)

## -----------------------------------------------------------------------------
select.samples.on.mutations(dataset.big, 100, desc = FALSE)

## -----------------------------------------------------------------------------
select.samples.on.mutations(dataset.big , 100, desc = FALSE) %>% select.genes.on.mutations(100)

## -----------------------------------------------------------------------------
corrplot.genes(dataset)

## -----------------------------------------------------------------------------
corrplot.samples(dataset)

## -----------------------------------------------------------------------------
# groups and counts equal genotypes
compactedDataset <- compact.dataset.easy(dataset)

## ---- echo=F------------------------------------------------------------------
head(compactedDataset)

## -----------------------------------------------------------------------------
samples <- as.matrix(compactedDataset %>% select(-freq))

## ----echo=F-------------------------------------------------------------------
samples

## -----------------------------------------------------------------------------
genes = colnames(samples)

## ----echo=F-------------------------------------------------------------------
genes

## -----------------------------------------------------------------------------
freqs = as.matrix(compactedDataset %>% ungroup() %>% select(freq))
freqs = freqs/sum(freqs)
freqs = c(freqs,0)

## ----echo=F-------------------------------------------------------------------
freqs

## -----------------------------------------------------------------------------
# prepare node labels listing the mutated genes for each node
labels <- prepare.labels(samples, genes)
# fix Colonal genotype absence, if needed
fix <- fix.clonal.genotype(samples, freqs, labels)
samples = fix[["samples"]]
freqs = fix[["freqs"]]
labels = fix[["labels"]]

## ----echo=F-------------------------------------------------------------------
samples

## -----------------------------------------------------------------------------
# compute edges based on subset relation
edges <- build.topology.subset(samples)

## -----------------------------------------------------------------------------
# remove transitive edges and prepare igraph object
g <- build.subset.graph(edges, labels)

## ---- echo=F------------------------------------------------------------------
V(g)$vertex.size <- rep(10, length(V(g)))
plot(g, vertex.size=rep(55, length(V(g))))

## -----------------------------------------------------------------------------
A <- as.matrix(as_adj(g))

## ---- echo=F------------------------------------------------------------------
A

## -----------------------------------------------------------------------------
no.of.children <- get.no.of.children(A,g)

## ---- echo=F------------------------------------------------------------------
no.of.children

## -----------------------------------------------------------------------------
upWeights <- computeUPW(g, freqs, no.of.children, A)

## ---- echo=F------------------------------------------------------------------
upWeights

## -----------------------------------------------------------------------------
normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)

## ---- echo=F------------------------------------------------------------------
normUpWeights

## -----------------------------------------------------------------------------
downWeights <- computeDWNW(g, freqs, no.of.children, A, normUpWeights)

## ---- echo=F------------------------------------------------------------------
downWeights

## -----------------------------------------------------------------------------
normDownWeights <- normalizeDWNW(g, freqs, no.of.children, A, downWeights)

## ---- echo=F------------------------------------------------------------------
normDownWeights

## -----------------------------------------------------------------------------
draw.ggraph(g, normDownWeights, labels)

## -----------------------------------------------------------------------------
draw.networkD3(g, normDownWeights, labels)

## -----------------------------------------------------------------------------
draw.visNetwork(g, normDownWeights, labels)

## ---- echo=FALSE--------------------------------------------------------------
# run ALL

