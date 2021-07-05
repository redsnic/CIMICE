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
#' @param original_dataset if not NULL, recover matching samples from the rownames of the original dataset
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
draw_visNetwork <- function(g, W, labels, original_dataset = NULL){

    V(g)$label <- labels
    E(g)$label <- W

    gr <- igraph_to_networkD3(g)
    gr$nodes$label = unlist(labels, use.names=FALSE)

    if(!is.null(original_dataset)){
        matching_samples <- get_samples_lables_matches(gr, original_dataset)
    }

    if(is.null(original_dataset)){
        matching_samples <- ""
    }else{
        matching_samples <- paste0("<br></p></font><p>Matching:<br>", matching_samples)
    }
   
    nodes <- data.frame(
        id = gr$nodes$name,
        label = seq(1,length(gr$nodes$label)), #gr$nodes$label,
        shadow = TRUE,
        shape = "box",
        title = paste0("<font color=\"black\"><p>Node:<br>",
                    gr$nodes$label,
                    matching_samples)
    )
    # compute genotype difference from source to destination
    diff <- map2(gr$links$source, gr$links$target,
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


#' Get associations between graph nodes and the samples of the original dataset
#' 
#' @param gr a networkD3 graph
#' @param original_dataset a dataset with row and column names (samples, genes)
#' 
#' @return a character vector of comma separated sample identifiers. 
#' These lists consists of the identifiers of the samples that match a certain node in the graph.
#' (follows node's orderd in the graph).
#'  
#' @examples
#' require(dplyr)
#' require(networkD3)
#' preproc <- example_dataset() %>% dataset_preprocessing
#' samples <- preproc[["samples"]]
#' freqs   <- preproc[["freqs"]]
#' labels  <- preproc[["labels"]]
#' genes   <- preproc[["genes"]]
#' g <- graph_non_transitive_subset_topology(samples, labels)
#' W <- compute_weights_default(g, freqs)
#' gr <- igraph_to_networkD3(g)
#' get_samples_lables_matches(gr, example_dataset())
#' 
get_samples_lables_matches <- function(gr, original_dataset){
    # prepare a matrix that will report for each pair (graph_node, sample) if the sample belongs to that graph_node 
    match_matrix <- matrix(rep(0, times=length(gr$nodes$name)*nrow(original_dataset)), ncol = nrow(original_dataset))
    # extract the list of considered genes in this analysis
    considered_genes <- strsplit( paste(gr$nodes$label, collapse = ", "), ",\\s*")[[1]] %>% unique 
    considered_genes <- considered_genes[considered_genes!="Clonal"]
    # find pairs (node, sample) with complete matching on wt/mutated genes
    for(node_i in seq(1,length(gr$nodes$name))){
        for(sample in seq(1,nrow(original_dataset))){
            matching <- TRUE
            present_genes <- strsplit(gr$nodes$label[[node_i]] , ",\\s*")[[1]]
            for(gene in considered_genes){
                if(original_dataset[sample,gene] != gene %in% present_genes){
                    matching <- FALSE
                }
            }
            match_matrix[node_i, sample] <- matching
        }
    }
    # prepare strings
    matching_samples <- rep("", times = length(gr$nodes$name))
    for(node_i in seq(1,length(gr$nodes$name))){
        matching_samples[node_i] <- rownames(original_dataset)[which(match_matrix[node_i, ] == TRUE)] %>% paste(collapse = ", ")
    }
    matching_samples
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

