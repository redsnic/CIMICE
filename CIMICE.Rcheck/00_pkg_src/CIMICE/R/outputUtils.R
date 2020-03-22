
# Draws the output graph using ggplot
# g <- graph to be drawn
# W <- weights on edges
# labels <- node labels
# digits <- precision for edges' weights
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

# Draws the output graph using networkD3
# g <- graph to be drawn
# W <- weights on edges
# labels <- node labels
draw.networkD3 <- function(g, W, labels){

    V(g)$label <- labels
    E(g)$label <- W

    gr <- igraph_to_networkD3(g)
    gr$nodes$label = unlist(labels, use.names=FALSE)
    # Create force directed network plot
    p<-forceNetwork(Links = gr$links, Nodes = gr$nodes,
                    Source = 'source', Target = 'target', Value = "value",
                    NodeID = 'label', Group = 'name', zoom=T, arrows = T, opacityNoHover = 1)

    p
}

# Draws the output graph using visNetwork
# g <- graph to be drawn
# W <- weights on edges
# labels <- node labels
draw.visNetwork <- function(g, W, labels){

    V(g)$label <- labels
    E(g)$label <- W

    gr <- igraph_to_networkD3(g)
    gr$nodes$label = unlist(labels, use.names=FALSE)

    nodes <- data.frame(
        id = gr$nodes$name,
        label = 1:length(gr$nodes$label), #gr$nodes$label,
        shadow = T,
        shape = "box",
        title = paste0("<font color=\"black\"><p>Node:<br>", gr$nodes$label,"</p></font>")
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
        shadow = T,
        label = signif(gr$links$value, digits = 4),
        title = paste0("<font color=\"black\"><p>Edge:<br> from:",
                       gr$nodes$label[gr$links$source+1] ,"<br>to : ",
                       gr$nodes$label[gr$links$target+1] ,"<br>New mutations:<br>",
                       diff, "</p></font>")
    )

    grp <- toVisNetworkData(g, idToLabel = TRUE)

    visNetwork(nodes, edges)
}
