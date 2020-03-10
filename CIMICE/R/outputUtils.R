
draw.ggraph <- function(g, W, labels, digits = 4){

    V(g)$label <- labels
    print(labels)
    print(V(g))
    print(E(g))
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

    edges <- data.frame(
        from = gr$links$source+1,
        to = gr$links$target+1,
        arrows = "to",
        shadow = T,
        label = signif(gr$links$value, digits = 4),
        # Sarebbero da indicare i geni che variano da A a B se A->B
        title = paste0("<font color=\"black\"><p>Edge:<br> from:", gr$links$source+1 ,"<br>to : ", gr$links$target+1,"</p></font>")
    )

    grp <- toVisNetworkData(g, idToLabel = TRUE)

    visNetwork(nodes, edges)
}
