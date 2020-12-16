# CPMC

require(igraph)
require(purrr)
require(ggraph)
require(tidygraph)


cpmc <- setRefClass("cpmc",
            fields = list(P = "matrix", labels = "list"),
            methods = list(
                get.id = function(label){
                    find.label(label, labels)
                },
                neighbours = function(id){
                    which( P[id,] > 0 )
                },
                graph = function(){
                    g <- graph_from_adjacency_matrix(P, weighted = T, mode = "directed")
                    g <- set.vertex.attribute(g, "label", index=V(g), map_chr(labels,~paste(., collapse = ", ")) )
                    g
                },
                set.edge = function(source , destination, p, normalize=TRUE){
                    if(normalize){
                        remainings <- sum(P[source,]) - p
                        for(j in 1:ncol(P))
                            P[source, j] <<- P[source, j]*remainings
                    }
                    P[source,destination] <<- p
                },
                add.node = function(label, unique = TRUE){
                    if(unique){
                        if(is.not.unique.label(label,labels)){
                            return(FALSE)
                        }
                    }
                    P <<- cbind(P, rep(0, nrow(P)))
                    P <<- rbind(P, rep(0, ncol(P)))
                    labels <<- c(labels, label)
                    return(TRUE)
                },
                remove.node = function(id = 0, label = NULL){
                    if(id != 0){
                        # if index is out of range does nothing
                        P <<- matrix(P[-id,-id], ncol=ncol(P)-1)
                        labels[[id]] <<- NULL
                        TRUE
                    }else if(!is.null(label)){
                        pos <- find.label(label,labels)
                        if(pos != 0){
                            P <<- matrix(P[-id,-id], ncol=ncol(P)-1)
                            labels[[pos]] <<- NULL
                            TRUE
                        }else{
                            FALSE
                        }
                    }else{
                        FALSE
                    }
                },
                plot = print_DTMC <- function(){
                    g <- graph_from_adjacency_matrix(P, weighted = TRUE) %>%
                        as_tbl_graph()
                    labs <- map_chr(labels, ~ paste(., collapse = ", "))
                    g$labels <- labs

                    ggraph(g) + theme_bw() +
                        geom_node_point(aes(color=as.factor(labs)), size = 10) +
                        geom_edge_loop(aes(alpha = log(E(g)$weight),
                                           label = round(E(g)$weight, 3)),
                                       end_cap = circle(10, 'mm'),
                                       start_cap = circle(10, 'mm'),
                                       angle_calc = 'along',
                                       label_dodge = unit(2.5, 'mm'),
                                       repel = T, color="red", arrow = arrow(length = unit(4, 'mm'))) +
                        geom_edge_link(aes(alpha = log(E(g)$weight),
                                           label = round(E(g)$weight, 3)),
                                       end_cap = circle(10, 'mm'),
                                       angle_calc = 'along',
                                       label_dodge = unit(2.5, 'mm'),
                                       start_cap = circle(10, 'mm'),
                                       repel = T, color="red", arrow = arrow(length = unit(4, 'mm'))) +
                        geom_node_label(aes(label = labs), repel = T) +
                        theme(legend.position = "none")
                },
                size = function(){ ncol(P) },
                check = function(){
                    if(ncol(P) != nrow(P)){
                        print("Non square probability matrix")
                        FALSE
                    }
                    if(ncol(P) != length(labels)){
                        print("Nodes and labels number does not match")
                        FALSE
                    }
                    ok <- TRUE
                    for(i in 1:nrow(P)){
                        if(sum(P[i,]) != 1){
                            print(paste("Non stochastic probability matrix (line:", i, "sum:", sum(P[i,]),")"))
                            ok <- FALSE
                        }
                    }
                    ok
                }
            ))


char.array.equal <- function(a, b){
    if(length(a) != length(b))
        return(FALSE)
    a <- sort(a)
    b <- sort(b)
    for(i in 1:length(a)){
        if(a[i] != b[i]){
            return(FALSE)
        }
    }
    return(TRUE)
}

is.not.unique.label <- function(label, labels){
    any(map_lgl(labels ,
                ~ char.array.equal(., label)))
}

find.label <- function(label, labels){
    ans <- which(TRUE == map_lgl(labels ,
                ~ char.array.equal(., label)))
    if(length(ans) == 0){
        0
    }else{
        ans[1]
    }
}

cpmc.apply <- function(D, procedure){
    n <- ncol(D$P)
    was.visited <- rep(FALSE, n)
    i <- 1
    while(i <= length(was.visited)){
        if(!was.visited[i]){
            # i and was.visited updates are demanded to the procedure
            out <- procedure(list(D = D, was.visited = was.visited, current=i))
            was.visited <- out$was.visited
            i <- out$current
        }else{
            i <- i+1
        }
    }
}

all.to.sink <- function(D){
    D$add.node('sink_node')
    sink_id <- D$size()
    make.sink <- function(input){
        i <- input$current
        (input$D)$set.edge(i,sink_id,1)
        input$was.visited[i] <- TRUE
        list(was.visited = input$was.visited, current=i+1)
    }
    cpmc.apply(D, make.sink)
}

cpmc.remove.singleton <- function(D){
    delete.singleton <- function(input){
        i <- input$current
        if(sum(input$D$P[i,]) - input$D$P[i,i] <= 0 ){
            input$D$remove.node(i)
            input$was.visited <- input$was.visited[-i]
            list(was.visited = input$was.visited, current=i)
        }else{
            input$was.visited[i] <- TRUE
            list(was.visited = input$was.visited, current=i+1)
        }
    }
    cpmc.apply(D, delete.singleton)
}


cpmc.normalize <- function(D){
    normalize.node <- function(input){
        i <- input$current
        tot <- sum(input$D$P[i,])
        if(tot != 0){
            for(j in 1:input$D$size()){
                input$D$set.edge(i,j,input$D$P[i,j]/tot, normalize = FALSE)
            }
        }
        input$was.visited[i] <- TRUE
        list(was.visited = input$was.visited, current=i+1)
    }
    cpmc.apply(D, normalize.node)
}

cpmc.apply.on.visit <- function(D, procedure, id=0, label=NULL){
    if(id==0){
        id <- D$get.id(label)
    }
    queue <- c(id)
    n <- ncol(D$P)
    was.visited <- rep(FALSE, n)
    while(length(queue) > 0){
        current <- queue[length(queue)]
        queue <- queue[-length(queue)]
        out <- procedure(list(D = D, was.visited = was.visited, current=current))
        was.visited <- out$was.visited
        neig <- keep( D$neighbours(current) , ~ was.visited[.] == FALSE)
        queue <- c(queue, neig)
    }
    was.visited
}

cpmc.prune <- function(D, id=0, label=NULL){
    simple.visit <- function(input){
        input$was.visited[input$current] <- TRUE
        list(was.visited = input$was.visited)
    }
    out <- cpmc.apply.on.visit(D, simple.visit, id=id, label=label)
    rem <- which(map_lgl(out, ~ !.))
    walk(rem, ~ D$remove.node(id = .))
}

cpmc.apply.treatment <- function(D, treatment, treatments){

    # find specific treatment
    found <- NULL
    for(t in treatments){
        if(t$name == treatment){
            found <- t
            break
        }
    }
    if(is.null(found)){
        return(FALSE)
    }
    print(found)
    # check effects
    for(effect in found$rules){
        apply.treatment <- NULL
        if(effect$action == "kill"){
            apply.treatment <- function(effect){
                force(effect)
                function(input){
                    i <- input$current

                    if(D$get.id('sink_node') == 0){
                        D$add.node('sink_node')
                        sink <- D$get.id('sink_node')
                        D$set.edge(sink, sink, 1, normalize = F)
                        input$was.visited <- c(input$was.visited, TRUE)
                    }
                    sink <- D$get.id('sink_node')

                    # add arc to sink
                    if(effect$guard(D$labels[[i]])){
                        D$set.edge(i,sink,effect$eff)
                    }

                    input$was.visited[i] <- TRUE
                    list(was.visited = input$was.visited, current=i+1)
                }
            }
            apply.treatment <- apply.treatment(effect)

        }else if(effect$action == "mutate"){
            print(paste("ERROR: mode not implemented:", effect$action))
            return(FALSE)
        }else{
            print(paste("ERROR: mode not implemented:", effect$action))
            return(FALSE)
        }

        cpmc.apply(D, apply.treatment)
    }

    return(TRUE)
}




