#' Filter dataset by genes' mutation count
#'
#' Dataset filtering on genes, based on their mutation count
#'
#' @param mutmatrix input dataset (mutational matrix) to be reduced
#' @param n number of genes to be kept
#' @param desc TRUE: select the n least mutated genes,
#' FALSE: select the n most mutated genes
#'
#' @return the modified dataset (mutational matrix)
#'
#' @examples
#'
#' # keep information on the 100 most mutated genes
#' select_genes_on_mutations(example_dataset(), 5)
#' # keep information on the 100 least mutated genes
#' select_genes_on_mutations(example_dataset(), 5, desc = FALSE)
#'
#' @export select_genes_on_mutations
select_genes_on_mutations <- function(mutmatrix, n, desc = TRUE){
    # transpose dataset to operate on genes
    temp.mutmatrix <- t(mutmatrix)
    # sort by sum
    sums <- rowSums(temp.mutmatrix, na.rm = TRUE)
    
    # sort and keep first n elements
    temp.mutmatrix <- temp.mutmatrix[order(sums, decreasing = desc),]
    temp.mutmatrix <- temp.mutmatrix[seq(1,min(n,nrow(temp.mutmatrix))),]
    
    # transpose back
    t(temp.mutmatrix)
}

#' Filter dataset by samples' mutation count
#'
#' Dataset filtering on samples, based on their mutation count
#'
#' @param mutmatrix input dataset (mutational matrix) to be reduced
#' @param n number of samples to be kept
#' @param desc T: select the n least mutated samples,
#' F: select the n most mutated samples
#'
#' @return the modified dataset (mutational matrix)
#'
#' @examples
#' require(dplyr)
#' # keep information on the 5 most mutated samples
#' select_samples_on_mutations(example_dataset(), 5)
#' # keep information on the 5 least mutated samples
#' select_samples_on_mutations(example_dataset(), 5, desc = FALSE)
#' # combine selections
#' select_samples_on_mutations(example_dataset() , 5, desc = FALSE) %>%
#'     select_genes_on_mutations(5)
#'
#' @export select_samples_on_mutations
select_samples_on_mutations <- function(mutmatrix, n, desc = TRUE){
    # sort by sum
    sums <- rowSums(mutmatrix, na.rm = TRUE)
    
    # sort and keep first n elements
    mutmatrix <- mutmatrix[order(sums, decreasing = desc),]
    mutmatrix <- mutmatrix[seq(1,min(n,nrow(mutmatrix))),]

    # transpose back
    mutmatrix
}

#' Radix sort for a binary matrix
#' 
#' Sort the rows of a binary matrix in ascending order
#' 
#' @param mat a binary matrix (of 0 and 1)
#' 
#' @return the sorted matrix
#' 
#' @examples 
#' require(Matrix)
#' m <- Matrix(c(1,1,0,1,0,0,0,1,1), sparse = TRUE, ncol = 3)
#' binary_radix_sort(m)
#' 
#' @export binary_radix_sort
binary_radix_sort <- function(mat){
    zeroes <- c()
    ones <- c()
    for(j in seq(ncol(mat), 1)){
        for(i in seq(1,nrow(mat))){
            if(mat[i,j] == 0){
                zeroes <- c(zeroes, i)
            }else{
                ones <- c(ones, i)
            }
        }
        mat <- mat[c(zeroes, ones), ]
        zeroes <- c()
        ones <- c()
    }
    mat
}


#' Compact dataset rows 
#'
#' Count duplicate rows
#' and compact the dataset (mutational). The column
#' 'freq' will contain the counts for each row.
#'
#' @param mutmatrix input dataset (mutational matrix)
#'
#' @return the compacted dataset (mutational matrix)
#'
#' @examples
#' compact_dataset(example_dataset())
#'
#' @export compact_dataset
compact_dataset <- function(mutmatrix){
    mutmatrix <- binary_radix_sort(mutmatrix)
    
    cnt <- 1
    counts <- c()
    # keep only a line per genotype (the first line is always new)
    valid_indexes <- c(1)
    for(i in seq(1,nrow(mutmatrix)-1)){
        if( all(mutmatrix[i,] == mutmatrix[i+1,]) ){
            cnt <- cnt + 1
        }else{
            valid_indexes <- c(valid_indexes, i+1)
            counts <- c(counts, cnt)
            cnt <- 1
        }
    }
    counts <- c(counts, cnt) 
    
    list(matrix = mutmatrix[c(valid_indexes), ], counts = counts)
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
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' build_topology_subset(samples)
#'
#' @export build_topology_subset
build_topology_subset <- function(samples){
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

#' Prepare node labels based on genotypes
#'
#' Prepare node labels so that each node is labelled with a
#' comma separated list of the alterated genes
#' representing its associated genotype.
#'
#' Note that after this procedure the user is
#' expected also to run fix_clonal_genotype
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
#' 
#' # compact
#' compactedDataset <- compact_dataset(example_dataset())
#' samples <- compactedDataset$matrix
#' 
#' # save genes' names
#' genes <- colnames(compactedDataset$matrix)
#' 
#' # keep the information on frequencies for further analysis
#' freqs <- compactedDataset$counts/sum(compactedDataset$counts)
#' 
#' # prepare node labels listing the mutated genes for each node
#' labels <- prepare_labels(samples, genes)
#'
#' @export prepare_labels
prepare_labels <- function(samples, genes){
    # prepare labels with the alterated genes accordingly to node's bit vector
    labels = apply(samples, MARGIN = 1, FUN = function(x) genes[x==1] )
    # concatenate the genes' names
    labels = lapply(labels,
                    function (x){
                        if (length(x) > 0)
                            paste(x, collapse=", ")
                        else
                            "Clonal"
                    }
    )
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
#' 
#' # compact
#' compactedDataset <- compact_dataset(example_dataset())
#' samples <- compactedDataset$matrix
#' 
#' # save genes' names
#' genes <- colnames(compactedDataset$matrix)
#' 
#' # keep the information on frequencies for further analysis
#' freqs <- compactedDataset$counts/sum(compactedDataset$counts)
#' 
#' # prepare node labels listing the mutated genes for each node
#' labels <- prepare_labels(samples, genes)
#' 
#' # fix Colonal genotype absence, if needed
#' fix <- fix_clonal_genotype(samples, freqs, labels)
#' 
#' samples <- fix[["samples"]]
#' freqs <- fix[["freqs"]]
#' labels <- fix[["labels"]]
#' 
#' list("samples" = samples, "freqs" = freqs,
#'     "labels" = labels, "genes" = genes)
#'
#' @export fix_clonal_genotype
fix_clonal_genotype <- function(samples, freqs, labels){
    # if no clonal genotype is found
    if (!(0 %in% rowSums(samples))){
        # add a 0 frequency genotype without mutations to the mutational matrix
        samples = rbind(samples, map_dbl(seq(1,ncol(samples)), function(x) 0) )
        freqs = c(freqs,0)
        # update labels
        labels = c(labels,"Clonal")
    }
    list("labels" = labels, "samples" = samples, "freqs"=freqs)
}

#' Remove transitive edges and prepare graph
#'
#' Create a graph from the "build_topology_subset" edge list, so
#' that it respects the subset relation, omitting the transitive edges.
#'
#' @param edges edge list, built from "build_topology_subset"
#' @param labels list of node labels, to be paired with the graph
#'
#' @return a graph with the subset topology, omitting transitive edges
#'
#' @examples
#' require(dplyr)
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' edges <- build_topology_subset(samples)
#' g <- build_subset_graph(edges, labels)
#'
#' @export build_subset_graph
build_subset_graph <- function(edges, labels){
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
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' A <- as_adj(g)
#' get_no_of_children(A, g)
#'
#' @export get_no_of_children
get_no_of_children <- function(A,g){
    no.of.children <- numeric(length(V(g)))
    for( v in V(g) ){
        no.of.children[v] = sum(A[v,])
    }
    no.of.children
}


