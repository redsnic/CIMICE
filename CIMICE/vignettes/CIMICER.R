## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ---- echo=F------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- include = FALSE---------------------------------------------------------
show.df <- function(df, w, h){
    df[1:h,1:w]
}
size.df <- function(df){
    paste("ncol:", ncol(df), " - nrow:", nrow(df))
}

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
# Remove transitive edges on a graph
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

## ---- results = 'hide'--------------------------------------------------------
# Read input dataset in CAPRI/CAPRESE format
dataset.big <- read.CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))

## ---- echo=F------------------------------------------------------------------
dataset.big %>% show.df(6,6) 
dataset.big %>% size.df()

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

## ---- echo=FALSE--------------------------------------------------------------
dataset

## -----------------------------------------------------------------------------
gene.mutations.hist(dataset.big)

## -----------------------------------------------------------------------------
sample.mutations.hist(dataset.big, binwidth = 10)

## ---- results = 'hide'--------------------------------------------------------
select.genes.on.mutations(dataset.big, 100)

## ---- echo=FALSE--------------------------------------------------------------
temp <- select.genes.on.mutations(dataset.big, 100)
temp %>% show.df(6,6) 
temp %>% size.df()

## ---- results = 'hide'--------------------------------------------------------
select.samples.on.mutations(dataset.big, 100, desc = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
temp <- select.samples.on.mutations(dataset.big, 100, desc = FALSE)
temp %>% show.df(6,6) 
temp %>% size.df()

## ---- results = 'hide'--------------------------------------------------------
select.samples.on.mutations(dataset.big , 100, desc = FALSE) %>% select.genes.on.mutations(100)

## ---- echo=FALSE--------------------------------------------------------------
temp <- select.samples.on.mutations(dataset.big , 100, desc = FALSE) %>% select.genes.on.mutations(100)
temp %>% show.df(6,6) 
temp %>% size.df()

## -----------------------------------------------------------------------------
corrplot.genes(dataset)

## -----------------------------------------------------------------------------
corrplot.samples(dataset)

## -----------------------------------------------------------------------------
# groups and counts equal genotypes
compactedDataset <- compact.dataset.easy(dataset)

## ---- echo=F------------------------------------------------------------------
compactedDataset

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

## ---- echo=F, out.height="50%", out.width="50%"-------------------------------
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

## ---- results = 'hide'--------------------------------------------------------
draw.networkD3(g, normDownWeights, labels)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("networkD3.png")

## ---- results = 'hide'--------------------------------------------------------
draw.visNetwork(g, normDownWeights, labels)

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("visGraph.png")

## ---- echo=FALSE--------------------------------------------------------------
# run ALL

## -----------------------------------------------------------------------------
## R Under development (unstable) (2020-11-20 r79451)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 20.04.1 LTS
## 
## Matrix products: default
## BLAS:   /usr/local/lib/R/lib/libRblas.so
## LAPACK: /usr/local/lib/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=it_IT.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=it_IT.UTF-8        LC_COLLATE=it_IT.UTF-8    
##  [5] LC_MONETARY=it_IT.UTF-8    LC_MESSAGES=it_IT.UTF-8   
##  [7] LC_PAPER=it_IT.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=it_IT.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] ggraph_2.0.4     purrr_0.3.4      ggcorrplot_0.1.3 visNetwork_2.0.9
##  [5] networkD3_0.4    relations_0.6-9  igraph_1.2.6     tidyr_1.1.2     
##  [9] glue_1.4.2       ggplot2_3.3.2    dplyr_1.0.2      CIMICE_0.1.0    
## [13] BiocStyle_2.19.0
## 
## loaded via a namespace (and not attached):
##  [1] tidyselect_1.1.0    xfun_0.19           slam_0.1-47        
##  [4] sets_1.0-18         reshape2_1.4.4      graphlayouts_0.7.1 
##  [7] lattice_0.20-41     colorspace_2.0-0    vctrs_0.3.5        
## [10] generics_0.1.0      htmltools_0.5.0     viridisLite_0.3.0  
## [13] yaml_2.2.1          rlang_0.4.8         pillar_1.4.7       
## [16] withr_2.3.0         tweenr_1.0.1        lifecycle_0.2.0    
## [19] plyr_1.8.6          stringr_1.4.0       munsell_0.5.0      
## [22] gtable_0.3.0        htmlwidgets_1.5.2   evaluate_0.14      
## [25] labeling_0.4.2      knitr_1.30          tidygraph_1.2.0    
## [28] highr_0.8           Rcpp_1.0.5          scales_1.1.1       
## [31] BiocManager_1.30.10 magick_2.5.2        jsonlite_1.7.1     
## [34] farver_2.0.3        gridExtra_2.3       ggforce_0.3.2      
## [37] digest_0.6.27       stringi_1.5.3       bookdown_0.21      
## [40] ggrepel_0.8.2       polyclip_1.10-0     grid_4.1.0         
## [43] tools_4.1.0         magrittr_2.0.1      tibble_3.0.4       
## [46] cluster_2.1.0       crayon_1.3.4        pkgconfig_2.0.3    
## [49] Matrix_1.2-18       MASS_7.3-53         ellipsis_0.3.1     
## [52] rmarkdown_2.5       viridis_0.5.1       R6_2.5.0           
## [55] compiler_4.1.0

