# calcolo dei pesi "UP"
computeUPW <- function(g, freqs, no.of.children, A){
    # creazione della struttura dati per l'output
    upWeights <<- seq(-1, -1, length.out = length(E(g)))
    upWeights <- sapply( 1:length(E(g)), function(x) 
        computeUPW.aux(g, x, freqs, no.of.children, A) )
    upWeights
}

# calcolo dei pesi "UP" 
computeUPW.aux = function(g, edge, freqs, no.of.children, A) {
    # programmazione dinamica
    if (upWeights[edge] != -1) return(upWeights[edge])
    # prendo testa e coda dell'arco considerato
    source = as_ids(tail_of(g,edge))
    destination = as_ids(head_of(g,edge))

    D = no.of.children[source]
    # frequenza osservata della sorgente
    P = freqs[source]
    
    W <- 0
    # ricorsione sui predecessori
    if(length(which(A[,source] >= 1)) != 0){
        # considero gli archi
        edges.entering.in.source <- lapply( as.list(which(A[,source] >= 1)), 
                                            function (x) c(x,source) )
        # passo agli ID
        edges.entering.in.source <- get.edge.ids(g, unlist(edges.entering.in.source)) 
        # ricorsione
        W <- sum( sapply( edges.entering.in.source , 
                          function (x) computeUPW.aux(g, x, freqs, no.of.children, A) ) )
    }
    # calcolo della formula
    upWeights[edge] <<- (1/D) * (P+W)
    
    return(upWeights[edge])
}

# normalizzazione dei pesi UP
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


# calcolo dei pesi "DOWN"
computeDWNW <- function(g, freqs, no.of.children, A, normUpWeights){
    # creazione della struttura dati per l'output
    downWeights <<- seq(-1, -1, length.out = length(E(g)))
    downWeights <- sapply( 1:length(E(g)), function(x)
        computeDWNW.aux(g, x, freqs, no.of.children, A, normUpWeights) )
    downWeights
}

computeDWNW.aux = function(g, edge, freqs, no.of.children, A, normUpWeights) {
    # Programmazione dinamica
    if (downWeights[edge] > -1) return(downWeights[edge])
    
    source = as_ids(tail_of(g,edge))
    destination = as_ids(head_of(g,edge))
    
    UPW = normUpWeights[edge]
    P = freqs[destination]
    
    W <- 0
    # considera i casi in cui ho archi uscenti
    if(length(which(A[destination,] >= 1)) != 0){
        edges.exiting.destination <- lapply( as.list(which(A[destination,] >= 1)),
                                             function (x) c(destination,x) )
        edges.exiting.destination <- get.edge.ids(g, unlist(edges.exiting.destination)) 
        W <- sum( sapply( edges.exiting.destination , function (x)
            computeDWNW.aux(g, x, freqs, no.of.children, A, normUpWeights)))
    }
    # applica la formula per il calcolo
    downWeights[edge] <<- (UPW) * (P+W)
    return(downWeights[edge])
}

# normalizzazione dei pesi "DOWN"
normalizeDWNW <- function(g, freqs, no.of.children, A, downWeights){
    normDownWeights <- seq(-1, -1, length.out = length(E(g)))
    
    for ( v in V(g) ){
        if(length(which(A[v,] >= 1)) != 0){
            edges.exiting.source <- lapply( as.list(which(A[v,] >= 1)), function (x) c(v,x) )
            edges.exiting.source <- get.edge.ids(g, unlist(edges.exiting.source)) 
            normVal <- sum( sapply( edges.exiting.source , function(x) downWeights[x] ) )
            for(e in edges.exiting.source){
                if(normVal == 0){ # è l'unico arco
                    normDownWeights[e] <- 1
                }else{
                    normDownWeights[e] <- downWeights[e]/normVal
                }
            }
        }
    }
    
    normDownWeights
}
