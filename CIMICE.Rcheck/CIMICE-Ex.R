pkgname <- "CIMICE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "CIMICE-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CIMICE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("build.subset.graph")
### * build.subset.graph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: build.subset.graph
### Title: Remove transitive edges and prepare graph
### Aliases: build.subset.graph

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
edges <- build.topology.subset(samples)
g <- build.subset.graph(edges, labels)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("build.subset.graph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("build.topology.subset")
### * build.topology.subset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: build.topology.subset
### Title: Compute subset relation as edge list
### Aliases: build.topology.subset

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
build.topology.subset(samples)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("build.topology.subset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compact.dataset.easy")
### * compact.dataset.easy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compact.dataset.easy
### Title: Compact dataset rows with plyr
### Aliases: compact.dataset.easy

### ** Examples

compact.dataset.easy(example.dataset())




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compact.dataset.easy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("compute.weights.default")
### * compute.weights.default

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute.weights.default
### Title: Compute default weights
### Aliases: compute.weights.default

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
compute.weights.default(g, freqs)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute.weights.default", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeDWNW")
### * computeDWNW

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeDWNW
### Title: Down weights computation
### Aliases: computeDWNW

### ** Examples

require(dplyr)
require(igraph)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
# prepare adj matrix
A <- as.matrix(as_adj(g))
# pre-compute exiting edges from each node
no.of.children <- get.no.of.children(A,g)
upWeights <- computeUPW(g, freqs, no.of.children, A)
normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
computeDWNW(g, freqs, no.of.children, A, normUpWeights)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeDWNW", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("computeUPW")
### * computeUPW

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: computeUPW
### Title: Up weights computation
### Aliases: computeUPW

### ** Examples

require(dplyr)
require(igraph)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
# prepare adj matrix
A <- as.matrix(as_adj(g))
# pre-compute exiting edges from each node
no.of.children <- get.no.of.children(A,g)
computeUPW(g, freqs, no.of.children, A)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("computeUPW", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("corrplot.from.df")
### * corrplot.from.df

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: corrplot.from.df
### Title: Dataframe based correlation plot
### Aliases: corrplot.from.df

### ** Examples

corrplot.from.df(example.dataset())




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("corrplot.from.df", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("corrplot.genes")
### * corrplot.genes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: corrplot.genes
### Title: Gene based correlation plot
### Aliases: corrplot.genes

### ** Examples

corrplot.genes(example.dataset())




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("corrplot.genes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("corrplot.samples")
### * corrplot.samples

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: corrplot.samples
### Title: Sample based correlation plot
### Aliases: corrplot.samples

### ** Examples

corrplot.samples(example.dataset())




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("corrplot.samples", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dataset.preprocessing")
### * dataset.preprocessing

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dataset.preprocessing
### Title: Run CIMICE preprocessing
### Aliases: dataset.preprocessing

### ** Examples

require(dplyr)
example.dataset() %>% dataset.preprocessing




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dataset.preprocessing", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dataset.preprocessing.population")
### * dataset.preprocessing.population

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dataset.preprocessing.population
### Title: Run CIMICE preprocessing for poulation format dataset
### Aliases: dataset.preprocessing.population

### ** Examples

require(dplyr)
example.dataset.withFreqs() %>% dataset.preprocessing.population




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dataset.preprocessing.population", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("draw.ggraph")
### * draw.ggraph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: draw.ggraph
### Title: Ggplot graph output
### Aliases: draw.ggraph

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
W <- compute.weights.default(g, freqs)
draw.ggraph(g,W,labels,digit = 3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("draw.ggraph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("draw.networkD3")
### * draw.networkD3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: draw.networkD3
### Title: NetworkD3 graph output
### Aliases: draw.networkD3

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
W <- compute.weights.default(g, freqs)
draw.networkD3(g,W,labels)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("draw.networkD3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("draw.visNetwork")
### * draw.visNetwork

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: draw.visNetwork
### Title: VisNetwork graph output (default)
### Aliases: draw.visNetwork

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
W <- compute.weights.default(g, freqs)
draw.visNetwork(g,W,labels)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("draw.visNetwork", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("example.dataset")
### * example.dataset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: example.dataset
### Title: Creates a simple example dataset
### Aliases: example.dataset

### ** Examples

example.dataset()




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("example.dataset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("example.dataset.withFreqs")
### * example.dataset.withFreqs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: example.dataset.withFreqs
### Title: Creates a simple example dataset with frequency column
### Aliases: example.dataset.withFreqs

### ** Examples

example.dataset.withFreqs()




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("example.dataset.withFreqs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fix.clonal.genotype")
### * fix.clonal.genotype

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fix.clonal.genotype
### Title: Manage Clonal genotype in data
### Aliases: fix.clonal.genotype

### ** Examples

require(dplyr)
# compact
compactedDataset <- compact.dataset.easy(example.dataset())
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
# fix clonal
fix <- fix.clonal.genotype(samples, freqs, labels)
samples = fix[["samples"]]
freqs = fix[["freqs"]]
labels = fix[["labels"]]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fix.clonal.genotype", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gene.mutations.hist")
### * gene.mutations.hist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gene.mutations.hist
### Title: Histogram of genes' frequencies
### Aliases: gene.mutations.hist

### ** Examples

gene.mutations.hist(example.dataset(), binwidth = 10)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gene.mutations.hist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get.no.of.children")
### * get.no.of.children

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get.no.of.children
### Title: Get number of children
### Aliases: get.no.of.children

### ** Examples

require(dplyr)
require(igraph)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
A <- as.matrix(as_adj(g))
get.no.of.children(A, g)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get.no.of.children", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("graph.non.transitive.subset.topology")
### * graph.non.transitive.subset.topology

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: graph.non.transitive.subset.topology
### Title: Default preparation of graph topology
### Aliases: graph.non.transitive.subset.topology

### ** Examples

require(dplyr)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
graph.non.transitive.subset.topology(samples, labels)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("graph.non.transitive.subset.topology", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("make.dataset")
### * make.dataset

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: make.dataset
### Title: Dataset line by line construction: initialization
### Aliases: make.dataset

### ** Examples

make.dataset(APC,P53,KRAS)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("make.dataset", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("normalizeDWNW")
### * normalizeDWNW

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: normalizeDWNW
### Title: Down weights normalization
### Aliases: normalizeDWNW

### ** Examples

require(dplyr)
require(igraph)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
# prepare adj matrix
A <- as.matrix(as_adj(g))
# pre-compute exiting edges from each node
no.of.children <- get.no.of.children(A,g)
upWeights <- computeUPW(g, freqs, no.of.children, A)
normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
downWeights <- computeDWNW(g, freqs, no.of.children, A, normUpWeights)
normalizeUPW(g, freqs, no.of.children, A, downWeights)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("normalizeDWNW", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("normalizeUPW")
### * normalizeUPW

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: normalizeUPW
### Title: Up weights normalization
### Aliases: normalizeUPW

### ** Examples

require(dplyr)
require(igraph)
preproc <- example.dataset() %>% dataset.preprocessing
samples <- preproc[["samples"]]
freqs   <- preproc[["freqs"]]
labels  <- preproc[["labels"]]
genes   <- preproc[["genes"]]
g <- graph.non.transitive.subset.topology(samples, labels)
# prepare adj matrix
A <- as.matrix(as_adj(g))
# pre-compute exiting edges from each node
no.of.children <- get.no.of.children(A,g)
upWeights <- computeUPW(g, freqs, no.of.children, A)
normalizeUPW(g, freqs, no.of.children, A, upWeights)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("normalizeUPW", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("prepare.labels")
### * prepare.labels

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: prepare.labels
### Title: Prepare node labels based on genotypes
### Aliases: prepare.labels

### ** Examples

require(dplyr)
# compact
compactedDataset <- compact.dataset.easy(example.dataset())
# remove frequencies column from the compacted dataset
samples <- as.matrix(compactedDataset %>% select(-freq))
# save genes' names
genes = colnames(samples)
# keep the information on frequencies for further analysis
freqs = as.matrix(compactedDataset %>% ungroup() %>% select(freq))
freqs = freqs/sum(freqs)
freqs = c(freqs,0)
# prepare node labels listing the mutated genes for each node
prepare.labels(samples, genes)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("prepare.labels", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("quick.run")
### * quick.run

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: quick.run
### Title: Run CIMICE defaults
### Aliases: quick.run

### ** Examples

quick.run(example.dataset())




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("quick.run", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read.CAPRI")
### * read.CAPRI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read.CAPRI
### Title: Read a "CAPRI" file
### Aliases: read.CAPRI

### ** Examples

#          "pathToDataset/myDataset.CAPRI"
read.CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read.CAPRI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read.CAPRI.string")
### * read.CAPRI.string

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read.CAPRI.string
### Title: Read "CAPRI" file from a string
### Aliases: read.CAPRI.string

### ** Examples

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




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read.CAPRI.string", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("read")
### * read

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: read
### Title: Read a "CAPRI" file
### Aliases: read

### ** Examples

read(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("read", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sample.mutations.hist")
### * sample.mutations.hist

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sample.mutations.hist
### Title: Histogram of samples' frequencies
### Aliases: sample.mutations.hist

### ** Examples

sample.mutations.hist(example.dataset(), binwidth = 10)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sample.mutations.hist", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("select.genes.on.mutations")
### * select.genes.on.mutations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: select.genes.on.mutations
### Title: Filter dataset by genes' mutation count
### Aliases: select.genes.on.mutations

### ** Examples


# keep information on the 100 most mutated genes
select.genes.on.mutations(example.dataset(), 5)
# keep information on the 100 least mutated genes
select.genes.on.mutations(example.dataset(), 5, desc = FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("select.genes.on.mutations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("select.samples.on.mutations")
### * select.samples.on.mutations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: select.samples.on.mutations
### Title: Filter dataset by samples' mutation count
### Aliases: select.samples.on.mutations

### ** Examples

require(dplyr)
# keep information on the 5 most mutated samples
select.samples.on.mutations(example.dataset(), 5)
# keep information on the 5 least mutated samples
select.samples.on.mutations(example.dataset(), 5, desc = FALSE)
# combine selections
select.samples.on.mutations(example.dataset() , 5, desc = FALSE) %>%
    select.genes.on.mutations(5)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("select.samples.on.mutations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("update_df")
### * update_df

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: update_df
### Title: Dataset line by line construction: add a sample
### Aliases: update_df

### ** Examples


require(dplyr)
make.dataset(APC,P53,KRAS)   %>%
    update_df("S1", 1, 0, 1) %>%
    update_df("S2", 1, 1, 1)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("update_df", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
