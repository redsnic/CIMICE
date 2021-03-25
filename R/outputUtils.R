utils::globalVariables(c("node1.label", "node2.label"))
#' ggplot graph output
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
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' W <- compute_weights_default(g, freqs)
#' draw_ggraph(g,W,labels,digit = 3)
#'
#' @export draw_ggraph
draw_ggraph <- function(g, W, labels, digits = 4){
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
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' W <- compute_weights_default(g, freqs)
#' draw_networkD3(g,W,labels)
#'
#' @export draw_networkD3
draw_networkD3 <- function(g, W, labels){
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
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' W <- compute_weights_default(g, freqs)
#' draw_visNetwork(g,W,labels)
#'
#' @export draw_visNetwork
draw_visNetwork <- function(g, W, labels){

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
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' W <- compute_weights_default(g, freqs)
#' to_dot(g,W,labels)
#'
#' @export to_dot
to_dot <- function(g, W, labels){

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

