# PD global structures
weight.env <- new.env()
weight.env$upWeights <- NULL
weight.env$downWeights <- NULL

#' Up weights computation
#'
#' Computes the up weights formula using
#' a Dinamic Programming approach (starting call),
#' see vignettes for further explaination.
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#'
#' @return a vector containing the Up weights for each edge
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
#' # prepare adj matrix
#' A <- as.matrix(as_adj(g))
#' # pre-compute exiting edges from each node
#' no.of.children <- get_no_of_children(A,g)
#' computeUPW(g, freqs, no.of.children, A)
#'
#' @export computeUPW
computeUPW <- function(g, freqs, no.of.children, A){
    # create data structure for Dynamic Programming
    weight.env$upWeights <- seq(-1, -1, length.out = length(E(g)))
    weight.env$upWeights <- map_dbl( seq(1,length(E(g))), function(x)
        computeUPW_aux(g, x, freqs, no.of.children, A) )
    weight.env$upWeights
}

#' Up weights computation (aux)
#'
#' Computes the up weights formula using
#' a Dinamic Programming approach (recursion),
#' see vignettes for further explaination.
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param edge the currently considered edge
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#'
#' @return a vector containing the Up weights for each edge
#'
computeUPW_aux = function(g, edge, freqs, no.of.children, A) {
    # check if this recursion was already computed
    if (weight.env$upWeights[edge] != -1) return(weight.env$upWeights[edge])
    # get source and destination of the currently considered edge
    source = as_ids(tail_of(g,edge))
    destination = as_ids(head_of(g,edge))
    # extract number of children of the source
    D = no.of.children[source]
    # observed genotype frequency of the source
    P = freqs[source]

    # compute weight
    W <- 0
    # recursion on predecessor
    if(length(which(A[,source] >= 1)) != 0){
        # consider edges entering in the source node
        edges.entering.in.source <-
            lapply( as.list(which(A[,source] >= 1)),
                    function (x) c(x,source) )
        # get edge IDs
        edges.entering.in.source <-
            get.edge.ids(g, unlist(edges.entering.in.source))
        # recurr
        W <- sum(map_dbl(edges.entering.in.source ,
                        function (x)
                            computeUPW_aux(g, x, freqs, no.of.children, A)))
    }
    # compute formula
    weight.env$upWeights[edge] <- (1/D) * (P+W)

    return(weight.env$upWeights[edge])
}

#' Up weights normalization
#'
#' Normalizes up weights so that the sum
#' of weights of edges entering in a node is 1
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#' @param upWeights Up weights as computed by computeUPW
#'
#' @return a vector containing the normalized Up weights for each edge
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
#' # prepare adj matrix
#' A <- as.matrix(as_adj(g))
#' # pre-compute exiting edges from each node
#' no.of.children <- get_no_of_children(A,g)
#' upWeights <- computeUPW(g, freqs, no.of.children, A)
#' normalizeUPW(g, freqs, no.of.children, A, upWeights)
#'
#' @export normalizeUPW
normalizeUPW <- function(g, freqs, no.of.children, A, upWeights) {
    normUpWeights <- seq(-1, -1, length.out = length(E(g)))

    for ( v in V(g) ){
        # almeno un arco entrante
        if(length(which(A[,v] >= 1)) != 0){
            edges.entering.in.source <-
                lapply( as.list(which(A[,v] >= 1)), function (x) c(x,v) )
            edges.entering.in.source <-
                get.edge.ids(g, unlist(edges.entering.in.source))
            normVal <-
                sum( map_dbl( edges.entering.in.source ,
                            function (x) upWeights[x]))
            for(e in edges.entering.in.source){
                if(normVal == 0){ # è l'unico arco
                    normUpWeights[e] <- 1
                }else{
                    normUpWeights[e] <- upWeights[e]/normVal
                }
            }
        }
    }
    normUpWeights
}

#' Down weights computation
#'
#' Computes the Down weights formula using
#' a Dinamic Programming approach (starting call),
#' see vignettes for further explaination.
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#' @param normUpWeights normalized up weights as computed by normalizeUPW
#'
#' @return a vector containing the Up weights for each edge
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
#' # prepare adj matrix
#' A <- as.matrix(as_adj(g))
#' # pre-compute exiting edges from each node
#' no.of.children <- get_no_of_children(A,g)
#' upWeights <- computeUPW(g, freqs, no.of.children, A)
#' normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
#' computeDWNW(g, freqs, no.of.children, A, normUpWeights)
#'
#' @export computeDWNW
computeDWNW <- function(g, freqs, no.of.children, A, normUpWeights){
    # create Dynamic Programming data structure
    weight.env$downWeights <- seq(-1, -1, length.out = length(E(g)))
    weight.env$downWeights <- map_dbl( seq(1,length(E(g))) , function(x)
        computeDWNW_aux(g, x, freqs, no.of.children, A, normUpWeights) )
    weight.env$downWeights
}

#' Down weights computation (aux)
#'
#' Computes the Down weights formula using
#' a Dinamic Programming approach (recursion),
#' see vignettes for further explaination.
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param edge the currently considered edge
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#' @param normUpWeights normalized up weights as computed by normalizeUPW
#'
#' @return a vector containing the Up weights for each edge
#'
computeDWNW_aux = function(g, edge, freqs, no.of.children, A, normUpWeights) {
    # check if this recursion was already computed
    if (weight.env$downWeights[edge] > -1) return(weight.env$downWeights[edge])
    # get source and destination of the currently considered edge
    source = as_ids(tail_of(g,edge))
    destination = as_ids(head_of(g,edge))

    UPW = normUpWeights[edge]
    P = freqs[destination]

    W <- 0
    # consider only nodes with exiting edges
    if(length(which(A[destination,] >= 1)) != 0){
        edges.exiting.destination <-
            lapply( as.list(which(A[destination,] >= 1)),
                    function (x) c(destination,x) )
        edges.exiting.destination <-
            get.edge.ids(g, unlist(edges.exiting.destination))
        W <- sum(map_dbl(
            edges.exiting.destination ,
            function (x)
                computeDWNW_aux(g, x, freqs, no.of.children, A, normUpWeights)))
    }
    # compute formula
    weight.env$downWeights[edge] <- (UPW) * (P+W)
    return(weight.env$downWeights[edge])
}

#' Down weights normalization
#'
#' Normalizes Down weights so that the sum
#' of weights of edges exiting a node is 1
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#' @param downWeights Down weights as computed by computeDWNW
#'
#' @return a vector containing the normalized Down weights for each edge
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
#' # prepare adj matrix
#' A <- as.matrix(as_adj(g))
#' # pre-compute exiting edges from each node
#' no.of.children <- get_no_of_children(A,g)
#' upWeights <- computeUPW(g, freqs, no.of.children, A)
#' normUpWeights <- normalizeUPW(g, freqs, no.of.children, A, upWeights)
#' downWeights <- computeDWNW(g, freqs, no.of.children, A, normUpWeights)
#' normalizeUPW(g, freqs, no.of.children, A, downWeights)
#'
#' @export normalizeDWNW
normalizeDWNW <- function(g, freqs, no.of.children, A, downWeights){
    normDownWeights <- seq(-1, -1, length.out = length(E(g)))

    for ( v in V(g) ){
        if(length(which(A[v,] >= 1)) != 0){
            edges.exiting.source <-
                lapply( as.list(which(A[v,] >= 1)), function (x) c(v,x))
            edges.exiting.source <-
                get.edge.ids(g, unlist(edges.exiting.source))
            normVal <-
                sum(map_dbl(edges.exiting.source,
                            function(x) downWeights[x]))
            for(e in edges.exiting.source){
                if(normVal == 0){
                    # if this is the only edge
                    normDownWeights[e] <- 1
                }else{
                    normDownWeights[e] <- downWeights[e]/normVal
                }
            }
        }
    }
    normDownWeights
}
