utils::globalVariables(c("sums","Samples"))

#' Read "CAPRI" file from a string
#'
#' Read a "CAPRI" formatted file, from a text string
#'
#' @param txt string in valid "CAPRI" format
#'
#' @return the described mutational matrix as a dataframe
#'
#' @examples
#' read.CAPRI.string("
#' s\\g A B C D
#' S1 0 0 0 1
#' S2 1 0 0 0
#' S3 1 0 0 0
#' S4 1 0 0 1
#' S5 1 1 0 1
#' S6 1 1 0 1
#' S7 1 0 1 1
#' S8 1 1 0 1
#' ")
#'
#' @export read.CAPRI.string
read.CAPRI.string <- function(txt){
    df <- read.csv(text=txt, sep = " ", row.names=NULL)
    colnames(df)[1] = "Samples"
    # manage repeated rownames
    tryCatch({
        rownames(df) = df[,1]
    },
    error = function(e){
        rownames(df) = seq(1,length(df[,1]))
    })
    df <- df %>% select(-Samples)
    df
}

#' Read a "CAPRI" file
#'
#' Read a "CAPRI" formatted file from the file system
#'
#' @param filename path to file
#'
#' @return the described mutational matrix as a dataframe
#'
#' @examples
#' #          "pathToDataset/myDataset.CAPRI"
#' read.CAPRI(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))
#'
#' @export read.CAPRI
read.CAPRI <- function(filename){
    df <- read.csv(filename, sep = " ", row.names=NULL)
    colnames(df)[1] = "Samples"
    # manage repeated rownames
    tryCatch({
        rownames(df) = df[,1]
    },
    error = function(e){
        rownames(df) = seq(1,length(df[,1]))
    })
    df <- df %>% select(-Samples)
    df
}

#' Read a "CAPRI" file
#'
#' Read a "CAPRI" formatted file, as read.CAPRI
#'
#' @param filename path to file
#'
#' @return the described mutational matrix as a dataframe
#'
#' @examples
#' read(system.file("extdata", "example.CAPRI", package = "CIMICE", mustWork = TRUE))
#'
#' @export read
read <- read.CAPRI

#' Dataset line by line construction: initialization
#'
#' Initialize a dataset for "line by line" creation
#'
#' @param ... gene names (do not use '"', the input
#' is automatically converted to strings)
#'
#' @return a mutational matrix without samples structured as dataframe
#'
#' @examples
#' make.dataset(APC,P53,KRAS)
#'
#' @export make.dataset
make.dataset <- function(...){
    df <- data.frame(matrix(ncol = length(enexprs(...)), nrow = 0))
    colnames(df) <- enexprs(...)
    df
}

#' Dataset line by line construction: add a sample
#'
#' Add a sample (a row) to an existing dataset.
#' This procedure is meant to be used with the "%>%" operator
#'
#' @param d an existing dataframe (mutational matrix)
#' @param sampleName the row (sample) name
#' @param ... sample's genotype (0/1 numbers)
#'
#' @return the modified dataframe (mutational matrix)
#'
#' @examples
#'
#' require(dplyr)
#' make.dataset(APC,P53,KRAS)   %>%
#'     update_df("S1", 1, 0, 1) %>%
#'     update_df("S2", 1, 1, 1)
#'
#'
#' @export update_df
update_df <- function(d, sampleName, ...){
    rnames <- rownames(d)
    cnames <- colnames(d)
    df <- rbind(d, unlist(enexprs(...)))
    rownames(df) <- c(rnames, sampleName)
    colnames(df) <- cnames
    df
}

#' Histogram of genes' frequencies
#'
#' Create the histogram of the genes' mutational frequencies
#'
#' @param df input dataframe (mutational matrix)
#' @param binwidth binwidth parameter for the histogram (as in ggplot)
#'
#' @examples
#' gene.mutations.hist(example.dataset(), binwidth = 10)
#'
#' @return the newly created histogram
#'
#' @export gene.mutations.hist
gene.mutations.hist <- function(df, binwidth = 1){
    # transpose dataframe
    temp.df <- as.data.frame(t(df[,seq(1,nrow(t(df)))]))
    # compute number of alteration per gene
    temp.df$sums <- rowSums(temp.df[,seq(1,ncol(temp.df))], na.rm = TRUE)
    # count number of genes
    nGenes <- nrow(temp.df)
    # create the histogram
    ggplot(temp.df, aes(x=sums)) +
        geom_histogram(color="black", fill="white", binwidth = binwidth) +
        labs(title = glue("Number of mutations per gene ({nGenes} genes)"),
            x = "Number of mutations", y = "Count")
}

#' Histogram of samples' frequencies
#'
#' Create the histogram of the samples' mutational frequencies
#'
#' @param df input dataframe (mutational matrix)
#' @param binwidth binwidth parameter for the histogram (as in ggplot)
#'
#' @examples
#' sample.mutations.hist(example.dataset(), binwidth = 10)
#'
#' @return the newly created histogram
#'
#' @export sample.mutations.hist
sample.mutations.hist <- function(df, binwidth = 1){
    temp.df <- df
    # compute number of alteration per sample
    temp.df$sums <- rowSums(temp.df[,seq(1,ncol(df))], na.rm = TRUE)
    # count number of samples
    nSample <- nrow(temp.df)
    # create the histogram
    ggplot(temp.df, aes(x=sums)) +
        geom_histogram(color="black", fill="white", binwidth = binwidth) +
        labs(title = glue("Number of mutations per sample ({nSample} samples)"),
            x = "Number of mutations",
            y = "Counts")
}

#' Filter dataset by genes' mutation count
#'
#' Dataset filtering on genes, based on their mutation count
#'
#' @param df input dataset (mutational matrix) to be reduced
#' @param n number of genes to be kept
#' @param desc TRUE: select the n least mutated genes,
#' FALSE: select the n most mutated genes
#'
#' @return the modified dataset (mutational matrix)
#'
#' @examples
#'
#' # keep information on the 100 most mutated genes
#' select.genes.on.mutations(example.dataset(), 5)
#' # keep information on the 100 least mutated genes
#' select.genes.on.mutations(example.dataset(), 5, desc = FALSE)
#'
#' @export select.genes.on.mutations
select.genes.on.mutations <- function(df, n, desc = TRUE){
    # transpose dataset to operate on genes
    temp.df <- as.data.frame(t(df))
    # sort by sum
    temp.df$sums <- rowSums(temp.df[,seq(1,ncol(temp.df))], na.rm = TRUE)
    if(desc){
        # decreasing sorting
        temp.df <- temp.df[order(-temp.df$sums),][seq(1,n),] %>% select(-sums)
    }else{
        # increasing sorting
        temp.df <- temp.df[order(temp.df$sums),][seq(1,n),] %>% select(-sums)
    }
    # transpose back and return to dataframe form
    temp.df <- as.data.frame(t(temp.df))
    # remove NA if the number of selected elements was too big
    temp.df %>% select_if(function(x) any(!is.na(x))) %>% drop_na()
}

#' Filter dataset by samples' mutation count
#'
#' Dataset filtering on samples, based on their mutation count
#'
#' @param df input dataset (mutational matrix) to be reduced
#' @param n number of samples to be kept
#' @param desc T: select the n least mutated samples,
#' F: select the n most mutated samples
#'
#' @return the modified dataset (mutational matrix)
#'
#' @examples
#' require(dplyr)
#' # keep information on the 5 most mutated samples
#' select.samples.on.mutations(example.dataset(), 5)
#' # keep information on the 5 least mutated samples
#' select.samples.on.mutations(example.dataset(), 5, desc = FALSE)
#' # combine selections
#' select.samples.on.mutations(example.dataset() , 5, desc = FALSE) %>%
#'     select.genes.on.mutations(5)
#'
#' @export select.samples.on.mutations
select.samples.on.mutations <- function(df, n, desc = TRUE){
    temp.df <- df
    # sort by sum
    temp.df$sums <- rowSums(temp.df[,2:ncol(temp.df)], na.rm = TRUE)
    if(desc){
        temp.df <- temp.df[order(-temp.df$sums),][seq(1,n),] %>% select(-sums)
    }else{
        temp.df <- temp.df[order(temp.df$sums),][seq(1,n),] %>% select(-sums)
    }
    # remove NA if the number of selected elements was too big
    temp.df %>% select_if(function(x) any(!is.na(x))) %>% drop_na()
}

#' Dataframe based correlation plot
#'
#' Prepare correlation plot based on a dataframe
#'
#' @param df input dataset
#'
#' @return the computed correlation plot
#'
#' @examples
#' corrplot.from.df(example.dataset())
#'
#' @export corrplot.from.df
corrplot.from.df <- function(df){
    # prepare correlation matrix
    corr <- round(cor(df), 1)
    corr[is.nan(corr)] <- 0
    corr[is.na(corr)] <- 0
    corr[corr == Inf] <- 0
    # prepare corrplot
    ggcorrplot(corr, hc.order = TRUE)
}

#' Gene based correlation plot
#'
#' Prepare correlation plot based on a
#' dataframe computed from genes' perspective
#'
#' @param df input dataset (mutational matrix)
#'
#' @return the computed correlation plot
#'
#' @examples
#' corrplot.genes(example.dataset())
#'
#' @export corrplot.genes
corrplot.genes <- corrplot.from.df

#' Sample based correlation plot
#'
#' Prepare correlation plot based on a dataframe
#' computed from samples' perspective
#'
#' @param df input dataset (mutational matrix)
#'
#' @return the computed correlation plot
#'
#' @examples
#' corrplot.samples(example.dataset())
#'
#' @export corrplot.samples
corrplot.samples <- function (df)
    corrplot.from.df(as.data.frame(t(df)))

#' Compact dataset rows with plyr
#'
#' Use plyr library to count duplicate rows
#' and compact the dataset (mutational). The column
#' 'freq' will contain the counts for each row.
#'
#' @param dataset input dataset (mutational matrix)
#'
#' @return the compacted dataset (mutational matrix)
#'
#' @examples
#' compact.dataset.easy(example.dataset())
#'
#' @export compact.dataset.easy
compact.dataset.easy <- function(dataset){
    dataset %>% plyr::count( colnames(dataset) )
}

#' Compute subset relation as edge list
#'
#' Create an edge list E representing the
#' 'subset' relation for binary strings so that:
#' \deqn{ (A,B) in E <=> forall(i) : A[i] -> B[i] }
#'
#' @param samples input dataset (mutational matrix) as matrix
#'
#' @return the computed edge list
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' build.topology.subset(samples)
#'
#' @export build.topology.subset
build.topology.subset <- function(samples){
    # computing subset relation
    edges = list()
    index = 1
    # simple for loop that computes
    # the subset relation pairwisely
    for(i in seq(1,nrow(samples))){
        for(j in seq(1,nrow(samples))){
            if(i!=j){
                r1 = samples[i,]
                r2 = samples[j,]
                if( ! ((-1) %in% (r1 - r2)) ) {
                    edges[[index]] <- c(j,i)
                    index <- index + 1
                }
            }
        }
    }
    edges
}

#' Remove transitive edges and prepare graph
#'
#' Create a graph from the "build.topology.subset" edge list, so
#' that it respects the subset relation, omitting the transitive edges.
#'
#' @param edges edge list, built from "build.topology.subset"
#' @param labels list of node labels, to be paired with the graph
#'
#' @return a graph with the subset topology, omitting transitive edges
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' edges <- build.topology.subset(samples)
#' g <- build.subset.graph(edges, labels)
#'
#' @export build.subset.graph
build.subset.graph <- function(edges, labels){
    # prepare the actual graph
    g = graph_from_edgelist(t(simplify2array(edges)))
    # remove transitive edges using transitive.reduction from the "nem" package
    # g = graph_from_adjacency_matrix(transitive.reduction(as.matrix(as_adj(g))))
    # altenative approach with "relations" package
    E <- transitive_reduction(
        endorelation(graph = as.list(data.frame(t(as_edgelist(g))))))
    g <- graph_from_adjacency_matrix(E$.Data$incidence)
    # add labels to node
    V(g)$label <- labels
    g
}

#' Prepare node labels based on genotypes
#'
#' Prepare node labels so that each node is labelled with a
#' comma separated list of the alterated genes
#' representing its associated genotype.
#'
#' Note that after this procedure the user is
#' expected also to run fix.clonal.genotype
#' to also add the clonal genortype to the
#' mutational matrix if it is not present.
#'
#' @param samples input dataset (mutational matrix) as matrix
#' @param genes list of gene names (in the columns' order)
#'
#' @return the computed edge list
#'
#' @examples
#' require(dplyr)
#' # compact
#' compactedDataset <- compact.dataset.easy(example.dataset())
#' # remove frequencies column from the compacted dataset
#' samples <- as.matrix(compactedDataset %>% select(-freq))
#' # save genes' names
#' genes = colnames(samples)
#' # keep the information on frequencies for further analysis
#' freqs = as.matrix(compactedDataset %>% ungroup() %>% select(freq))
#' freqs = freqs/sum(freqs)
#' freqs = c(freqs,0)
#' # prepare node labels listing the mutated genes for each node
#' prepare.labels(samples, genes)
#'
#' @export prepare.labels
prepare.labels <- function(samples, genes){
    # prepare labels with the alterated genes accordingly to node's bit vector
    labels = apply(samples, MARGIN = 1, FUN = function(x) genes[x==1] )
    # concatenate the genes' names
    labels = lapply(labels, function (x) paste(x, collapse=", "))
    labels
}

#' Manage Clonal genotype in data
#'
#' Fix the absence of the clonal genotype in the data (if needed)
#'
#' @param samples input dataset (mutational matrix) as matrix
#' @param freqs genotype frequencies (in the rows' order)
#' @param labels list of gene names (in the columns' order)
#'
#' @return a named list containing the fixed "samples", "freqs" and "labels"
#'
#' @examples
#' require(dplyr)
#' # compact
#' compactedDataset <- compact.dataset.easy(example.dataset())
#' # remove frequencies column from the compacted dataset
#' samples <- as.matrix(compactedDataset %>% select(-freq))
#' # save genes' names
#' genes = colnames(samples)
#' # keep the information on frequencies for further analysis
#' freqs = as.matrix(compactedDataset %>% ungroup() %>% select(freq))
#' freqs = freqs/sum(freqs)
#' freqs = c(freqs,0)
#' # prepare node labels listing the mutated genes for each node
#' labels <- prepare.labels(samples, genes)
#' # fix clonal
#' fix <- fix.clonal.genotype(samples, freqs, labels)
#' samples = fix[["samples"]]
#' freqs = fix[["freqs"]]
#' labels = fix[["labels"]]
#'
#' @export fix.clonal.genotype
fix.clonal.genotype <- function(samples, freqs, labels){
    # if no clonal genotype is found
    if (!(0 %in% apply(samples,MARGIN=1, sum))){
        # add a 0 frequency genotype without mutations to the mutational matrix
        samples = rbind(samples, map_dbl(seq(1,ncol(samples)), function(x) 0) )
        freqs = c(freqs,0)
        # update labels
        labels = c(labels,"Clonal")
    }
    list("labels" = labels, "samples" = samples, "freqs"=freqs)
}

#' Get number of children
#'
#' Compute number of children for each node given an adj matrix
#'
#' @param A Adjacency matrix of the graph g
#' @param g a graph
#'
#' @return a vector containing the number of children for each node in g
#'
#' @examples
#' require(dplyr)
#' require(igraph)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' A <- as.matrix(as_adj(g))
#' get.no.of.children(A, g)
#'
#' @export get.no.of.children
get.no.of.children <- function(A,g){
    no.of.children <- numeric(length(V(g)))
    for( v in V(g) ){
        no.of.children[v] = sum(A[v,])
    }
    no.of.children
}


