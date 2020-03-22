#' Up weights computation
#'
#' Computes the up weights formula:
#'
#' \deqn{ W_{up}(<a,b> in E) = (1/|L_a|)(P(a) + \Sigma{for: x in \Pi_a}{of: W_{up}(<x,a>))}}
#'
#' using a Dinamic Programming approach (starting call), see vignettes for further explaination.
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#'
#' @return a vector containing the Up weights for each edge
#'
#' @examples
#' computeUPW(g, freqs, no.of.children, A)
#'
#' @export
computeUPW <- function(g, freqs, no.of.children, A){
    # create data structure for Dynamic Programming
    upWeights <<- seq(-1, -1, length.out = length(E(g)))
    upWeights <- sapply( 1:length(E(g)), function(x)
        computeUPW.aux(g, x, freqs, no.of.children, A) )
    upWeights
}

#' Up weights computation (aux)
#'
#' Computes the up weights formula:
#'
#' \deqn{ W_{up}(<a,b> in E) = (1/|L_a|)(P(a) + \Sigma{for: x in \Pi_a}{of: W_{up}(<x,a>))}}
#'
#' using a Dinamic Programming approach (recursion), see vignettes for further explaination.
#'
#' @param g graph (a Directed Acyclic Graph)
#' @param edge the currently considered edge
#' @param freqs observed genotype frequencies
#' @param no.of.children number of children for each node
#' @param A adjacency matrix of G
#'
#' @return a vector containing the Up weights for each edge
#'
#' @examples
#' computeUPW.aux(g, edge, freqs, no.of.children, A)
#'
computeUPW.aux = function(g, edge, freqs, no.of.children, A) {
    # Check if this recursion was already computed
    if (upWeights[edge] != -1) return(upWeights[edge])
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
        edges.entering.in.source <- lapply( as.list(which(A[,source] >= 1)),
                                            function (x) c(x,source) )
        # get edge IDs
        edges.entering.in.source <- get.edge.ids(g, unlist(edges.entering.in.source))
        # recurr
        W <- sum( sapply( edges.entering.in.source ,
                          function (x) computeUPW.aux(g, x, freqs, no.of.children, A) ) )
    }
    # compute formula
    upWeights[edge] <<- (1/D) * (P+W)

    return(upWeights[edge])
}

#' Up weights normalization
#'
#' Normalizes up weights so that the sum of weights of edges entering in a node is 1
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
#' normalizeUPW(g, freqs, no.of.children, A, upWeights)
#'
#' @export
normalizeUPW <- function(g, freqs, no.of.children, A, upWeights) {
    normUpWeights <- seq(-1, -1, length.out = length(E(g)))

    for ( v in V(g) ){
        # almeno un arco entrante
        if(length(which(A[,v] >= 1)) != 0){
            edges.entering.in.source <- lapply( as.list(which(A[,v] >= 1)), function (x) c(x,v) )
            edges.entering.in.source <- get.edge.ids(g, unlist(edges.entering.in.source))
            normVal <- sum( sapply( edges.entering.in.source , function (x) upWeights[x] ) )
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
#' Computes the Down weights formula:
#'
#' W_{down}(<a,b>) = NormalizedW_{up}(<a,b>)(P(b) + \Sigma{for: x in L_b}{of: W_{down}(<b,x>)}
#'
#' using a Dinamic Programming approach (starting call), see vignettes for further explaination.
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
#' computeDWNW(g, freqs, no.of.children, A, normUpWeights)
#'
#' @export
computeDWNW <- function(g, freqs, no.of.children, A, normUpWeights){
    # create Dynamic Programming data structure
    downWeights <<- seq(-1, -1, length.out = length(E(g)))
    downWeights <- sapply( 1:length(E(g)), function(x)
        computeDWNW.aux(g, x, freqs, no.of.children, A, normUpWeights) )
    downWeights
}

#' Up weights computation (aux)
#'
#' Computes the Down weights formula:
#'
#' W_{down}(<a,b>) = NormalizedW_{up}(<a,b>)(P(b) + \Sigma{for: x in L_b}{of: W_{down}(<b,x>)}
#'
#' using a Dinamic Programming approach (recursion), see vignettes for further explaination.
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
#' @examples
#' computeUPW.aux(g, edge, freqs, no.of.children, A)
#'
computeDWNW.aux = function(g, edge, freqs, no.of.children, A, normUpWeights) {
    # Check if this recursion was already computed
    if (downWeights[edge] > -1) return(downWeights[edge])
    # get source and destination of the currently considered edge
    source = as_ids(tail_of(g,edge))
    destination = as_ids(head_of(g,edge))

    UPW = normUpWeights[edge]
    P = freqs[destination]

    W <- 0
    # consider only nodes with exiting edges
    if(length(which(A[destination,] >= 1)) != 0){
        edges.exiting.destination <- lapply( as.list(which(A[destination,] >= 1)),
                                             function (x) c(destination,x) )
        edges.exiting.destination <- get.edge.ids(g, unlist(edges.exiting.destination))
        W <- sum( sapply( edges.exiting.destination , function (x)
            computeDWNW.aux(g, x, freqs, no.of.children, A, normUpWeights)))
    }
    # compute formula
    downWeights[edge] <<- (UPW) * (P+W)
    return(downWeights[edge])
}

#' Down weights normalization
#'
#' Normalizes down weights so that the sum of weights of edges exiting a node is 1
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
#' normalizeUPW(g, freqs, no.of.children, A, downWeights)
#'
#' @export
normalizeDWNW <- function(g, freqs, no.of.children, A, downWeights){
    normDownWeights <- seq(-1, -1, length.out = length(E(g)))

    for ( v in V(g) ){
        if(length(which(A[v,] >= 1)) != 0){
            edges.exiting.source <- lapply( as.list(which(A[v,] >= 1)), function (x) c(v,x))
            edges.exiting.source <- get.edge.ids(g, unlist(edges.exiting.source))
            normVal <- sum( sapply( edges.exiting.source , function(x) downWeights[x] ) )
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
