utils::globalVariables(c("node1.label", "node2.label"))
#' Ggplot graph output
#'
#' Draws the output graph using ggplot
#'
#' @param g graph to be drawn
#' @param W weights on edges
#' @param labels node labels
#' @param digits precision for edges' weights
#'
#' @return ggraph object representing g as described
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' W <- compute.weights.default(g, freqs)
#' draw.ggraph(g,W,labels,digit = 3)
#'
#' @export draw.ggraph
draw.ggraph <- function(g, W, labels, digits = 4){
    V(g)$label <- labels

    weights <- signif(W, digits = 4)

    ggraph(g, layout = "kk") +
        geom_node_point() +
        geom_node_label(aes(label=labels)) +
        geom_edge_link(
            aes(label=weights, start_cap = label_rect(node1.label),
                end_cap = label_rect(node2.label)),
            arrow = arrow(length = unit(4, 'mm')),
            label_dodge = unit(2.5, 'mm'),
            angle_calc = 'along'
        )
}

#' NetworkD3 graph output
#'
#' Draws the output graph using networkD3
#'
#' @param g graph to be drawn
#' @param W weights on edges
#' @param labels node labels
#'
#' @return networkD3 object representing g as described
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' W <- compute.weights.default(g, freqs)
#' draw.networkD3(g,W,labels)
#'
#' @export draw.networkD3
draw.networkD3 <- function(g, W, labels){

    V(g)$label <- labels
    E(g)$label <- W

    gr <- igraph_to_networkD3(g)
    gr$nodes$label = unlist(labels, use.names=FALSE)
    # Create force directed network plot
    p<-forceNetwork(Links = gr$links, Nodes = gr$nodes,
                    Source = 'source', Target = 'target', Value = "value",
                    NodeID = 'label', Group = 'name',
                    zoom=TRUE, arrows = TRUE, opacityNoHover = 1)
    p
}

#' VisNetwork graph output (default)
#'
#' Draws the output graph using VisNetwork
#'
#' @param g graph to be drawn
#' @param W weights on edges
#' @param labels node labels
#'
#' @return visNetwork object representing g as described
#'
#' @examples
#' require(dplyr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' W <- compute.weights.default(g, freqs)
#' draw.visNetwork(g,W,labels)
#'
#' @export draw.visNetwork
draw.visNetwork <- function(g, W, labels){

    V(g)$label <- labels
    E(g)$label <- W

    gr <- igraph_to_networkD3(g)
    gr$nodes$label = unlist(labels, use.names=FALSE)

    nodes <- data.frame(
        id = gr$nodes$name,
        label = seq(1,length(gr$nodes$label)), #gr$nodes$label,
        shadow = TRUE,
        shape = "box",
        title = paste0("<font color=\"black\"><p>Node:<br>",
                    gr$nodes$label,"</p></font>")
    )
    # compute genotype difference from source to destination
    diff = map2(gr$links$source, gr$links$target,
                function(x, y){
                    a <- unlist(strsplit(gr$nodes$label[x+1], ", "))
                    b <- unlist(strsplit(gr$nodes$label[y+1], ", "))
                    d <- setdiff(b,a)
                    d
                })

    edges <- data.frame(
        from = gr$links$source+1,
        to = gr$links$target+1,
        arrows = "to",
        shadow = TRUE,
        label = signif(gr$links$value, digits = 4),
        title = paste0("<font color=\"black\"><p>Edge:<br> from:",
                    gr$nodes$label[gr$links$source+1] ,"<br>to : ",
                    gr$nodes$label[gr$links$target+1] ,
                    "<br>New mutations:<br>",
                    diff, "</p></font>")
    )

    grp <- toVisNetworkData(g, idToLabel = TRUE)

    visNetwork(nodes, edges)
}


#' Dot graph output
#'
#' Represents this graph in dot format (a textual output format)
#'
#' @param g graph to be drawn
#' @param W weights on edges
#' @param labels node labels
#'
#' @return a string representing the graph in dot format
#'
#' @examples
#' require(dplyr)
#' require(purrr)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' W <- compute.weights.default(g, freqs)
#' to.dot(g,W,labels)
#'
#' @export to.dot
to.dot <- function(g, W, labels){

    V(g)$label <- labels
    E(g)$label <- W

    str <- paste(
        "digraph G{",
        paste(
            imap_chr(V(g)$label, ~ paste( .y, '[label="', .x ,'"]', sep="") ),
            collapse = "\n"),
        paste(
            imap_chr(E(g), ~ paste(tail_of(g, E(g)[[.y]]), " -> ", head_of(g, E(g)[[.y]]), '[label="', round(W[.y], digits = 3) , '"]', sep="" )),
            collapse = "\n"),
        "}",
        sep="\n"
    )

    str
}

#' to prism module
#'
#' convert this graph to a prism model
#'
#' @param g graph to be drawn
#' @param W weights on edges
#' @param labels node labels
#'
#' @return a string representing the graph as a prism model
#'
#' @examples
#' require(dplyr)
#' require(purrr)
#' require(igraph)
#' preproc <- example.dataset() %>% dataset.preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph.non.transitive.subset.topology(samples, labels)
#' W <- compute.weights.default(g, freqs)
#' to.prism(g,W,labels)
#'
#' @export to.prism
to.prism <- function(g, W, labels){
    V(g)$label <- labels
    E(g)$label <- W
    g_ <- set.edge.attribute(g, "weight", index=E(g), W)
    W_ <- as_adjacency_matrix(g_, attr = "weight",sparse = F)
    W_ <- add.self.loops(W_)
    str <- paste(
        "dtmc", "module M", paste("\tx : [1..", length(V(g)) ,"] init 1;", sep=""),
        paste(format_transition_probabilities(W_), collapse = "\n"),
        "endmodule",
        sep="\n"
    )

    str
}

#' format transition probabilities
#'
#' convert a transition probability matrix into transitions written
#' in prism syntax
#'
#' @param W a transition probability matrix
#'
#' @return a character vector containing the transitions in prism syntax
#'
#' @examples
#' W <- matrix(c(0.3,0.3,0.4,0,0.5,0.5,1,0,0), nrow=3, byrow=TRUE)
#' format_transition_probabilities(W)
#'
#' @export format_transition_probabilities
format_transition_probabilities <- function(W){
    out <- rep("", ncol(W))
    for(i in 1:nrow(W)){
        first <- TRUE
        for(j in 1:ncol(W)){
            if(W[i,j] != 0){
                if(first){
                    out[i] <- paste("\t[] x=", i, " -> ", W[i,j], ":(x'=", j, ")", sep="")
                    first <- FALSE
                }else{
                    out[i] <- paste(out[i], " + ", W[i,j], ":(x'=", j, ")", sep="")
                }
            }
        }
        out[i] <- paste(out[i], ";", sep="")
    }
    out
}

#' add self loops on sink nodes
#'
#' given a transition probability matrix, this procedures adds
#' self loops to sink nodes (nodes without any exiting edge)
#'
#' @param W a transition probability matrix
#'
#' @return the updated probability matrix with self loops
#'
#' @examples
#' W = matrix(c(0,0,0,0,0.5,0.5,1,0,0), nrow=3, byrow=TRUE)
#' add.self.loops(W)
#'
#' @export add.self.loops
add.self.loops <- function(W){
    for(i in 1:nrow(W)){
        # check if there are no exiting edges
        if(sum(W[i,]) == 0){
            W[i,i] <- 1
        }
    }
    W
}
