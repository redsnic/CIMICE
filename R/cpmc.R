# CPMC

require(igraph)
require(purrr)
require(ggraph)
require(tidygraph)


#' Class to represent Cancer Progression Markov Chains CPMC
#'
#' This representation is also well fit to encode any type of
#' labeled Markov chains
#'
#' @name cpmc
#' @field P : a transition probability matrix
#' @field labels : list of character vector representing genotypes
#'                  (or a list of atomic proposition)
#' @export cpmc
cpmc <- setRefClass("cpmc",
            fields = list(P = "matrix", labels = "list"),
            methods = list(
                #' @description
                #' get node id (index in P) form label
                #' @param label a label (character vector)
                #' @return the id of a label if found, 0 else
                get.id = function(label){
                    find.label(label, labels)
                },
                #' @description
                #' get array of the ids of neighbors nodes
                #' @param id node id
                #' @return array of ids of neighbor nodes
                neighbours = function(id){
                    which( P[id,] > 0 )
                },
                #' @description
                #' creates an igraph object representing the current cpmc
                #' @return the computed graph
                graph = function(){
                    g <- graph_from_adjacency_matrix(P, weighted = T, mode = "directed")
                    g <- set.vertex.attribute(g, "label", index=V(g), map_chr(labels,~paste(., collapse = ", ")) )
                    g
                },
                #' @description
                #' set probability value of a transition from a node to another
                #' @param source source node id
                #' @param destination destination node id
                #' @param p transition probability
                #' @param normalize set TRUE if the remaining probabilities
                #'        of exiting from the source must be normalized
                set.edge = function(source, destination, p, normalize=TRUE){
                    if(p > 0){
                        if(normalize){
                            remainings <- max(c(0, 1 - p))
                            total <- sum(P[source,]) - P[source,destination]
                            if(total > 0){
                                for(j in 1:ncol(P))
                                    P[source, j] <<- (P[source, j]/total)*remainings
                            }
                        }
                        P[source,destination] <<- p
                    }
                },
                #' @description
                #' add a node to this cpmc
                #' @param label character vector describing genotype or
                #'        atomic properties
                #' @param unique if TRUE check before adding the node if it
                #'        has a unique label
                #' @return TRUE if the operation was successful
                add.node = function(label, unique = TRUE){
                    if(unique){
                        if(is.not.unique.label(label,labels)){
                            return(FALSE)
                        }
                    }
                    P <<- cbind(P, rep(0, nrow(P)))
                    P <<- rbind(P, rep(0, ncol(P)))
                    labels <<- c(labels, list(label) )
                    return(TRUE)
                },
                #' @description
                #' remove a node from its id or label
                #' @param id node id
                #' @param label node label
                #'        (removes the node with this label and lowest id)
                #' @return TRUE if the operation was successful
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
                #' @description
                #' plot this cpmc
                #' @return the computed plot
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
                #' @description
                #' retrieve number of nodes of this cpmc
                #' @return number of nodes of this cpmc
                size = function(){ ncol(P) },
                #' @description
                #' check if this is a valid cpmc
                #' @return TRUE if all tests are passed
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
                },
                #' @description
                #' format this cpmc as a PRISM model
                #' @return a PRISM formatted string representing this cpmc
                to.prism = function(name, initial=NULL, initial_state_description = 0){
                    if(!.self$check()){
                        return(NULL)
                    }
                    # get initial state automatically if not specified
                    if(is.null(initial)){
                        initial <- .self$get.id("Clonal")
                        if(initial == 0){
                            initial <- 1
                        }
                    }
                    # auto-compute if initial_state_description is not defined
                    if(length(initial_state_description) == 1 && initial_state_description[1] == 0){
                        initial_state_description <- find.occ(c("Clonal"), all_labels <- labels %>% unlist() %>% unique())
                    }else{
                        initial_state_description <- find.occ(initial_state_description, all_labels <- labels %>% unlist() %>% unique())
                    }

                    lable.table <- make.label.table(.self$labels)
                    str <- paste(
                        "dtmc", paste("module ", name), paste("\tx : [1..", .self$size() ,"] init ", initial, ";", sep=""),
                        paste(labels_to_prism_properties(lable.table, initial_state_description), collapse = "\n"),
                        paste(format_transition_probabilities_with_labels(.self$P, lable.table), collapse = "\n"),
                        "endmodule",
                        sep="\n"
                    )
                    str
                }
            ))


#' Prepare PRISM header
#'
#' Prepare PRISM header by defining state variable and lable variables.
#'
#' @param lable.table matrix of the characteristic functions of the labels
#' @param initial_state_description list of labels that are true in the initial state
#' TODO: the initial state label could be retrieved directly by the labels of the cpmc
#' @return translated value
labels_to_prism_properties <- function(lable.table, initial_state_description){
    out <- c()
    for(i in 1:ncol(lable.table)){
        out <- c(out, paste("\t_", colnames(lable.table)[i], " : bool init ", bool.prism.translate(initial_state_description[i]) ,";", sep=""))
    }
    out
}

#' Boolean to PRISM Boolean
#'
#' Convert Boolean values in PRISM language
#'
#' @param  b boolean value
#' @return translated value
bool.prism.translate <- function(b){
    if(b){
        "true"
    }else{
        "false"
    }
}

#' Format transitions in prism language
#'
#' Format transitions in the prism language. Also manages property labels.
#'
#' @param P transition probability matrix
#' @param lable.table matrix of the characteristic function of labels
#' @return transitions in prism language
#'
#' @examples
#' cimice.out <- quick.run(example.dataset())
#' DTMC <- to.cpmc(cimice.out$topology, cimice.out$weights, cimice.out$labels)
#' lt <- make.label.table(DTMC$labels)
#' format_next_state(1,lt)
#'
#' @export format_transition_probabilities_with_labels
format_transition_probabilities_with_labels <- function(P , lable.table){
    out <- rep("", ncol(P))
    for(i in 1:nrow(P)){
        first <- TRUE
        for(j in 1:ncol(P)){
            if(P[i,j] != 0){
                if(first){
                    out[i] <- paste("\t[] x=", i, " -> ", P[i,j], ":(x'=", j, ")", format_next_state(j,lable.table),sep="")
                    first <- FALSE
                }else{
                    out[i] <- paste(out[i], "\n\t\t+ ", P[i,j], ":(x'=", j, ")", format_next_state(j,lable.table), sep="")
                }
            }
        }
        out[i] <- paste(out[i], ";", sep="")
    }
    out
}

#' Format labels in prism language
#'
#' Format property labels in the prism language. Labels are represented
#' by Boolean values, in the contest of tumor philogenetics, they represent
#' genotypes. Labels are updated by the transitions.
#'
#' @param dest_id destination node to be considered
#' @param lable.table matrix of the characteristic function of labels
#' @return transition labels in prism language
#'
#' @examples
#' cimice.out <- quick.run(example.dataset())
#' DTMC <- to.cpmc(cimice.out$topology, cimice.out$weights, cimice.out$labels)
#' lt <- make.label.table(DTMC$labels)
#' format_next_state(1,lt)
#'
#' @export format_next_state
format_next_state <- function(dest_id, label.table){
    out <- c()
    for(j in 1:ncol(label.table)){
        out <- c(out, paste(" & (", "_", colnames(label.table)[j],"'=", bool.prism.translate(label.table[dest_id, j]) ,")",sep=""))
    }
    paste(out, collapse = "")
}

#' Create table of the labels
#'
#' Translate labels to a matrix of characteristic functions. The
#' characteristic functions follow the ascending order of sorted labels.
#'
#' @param labels set of labels
#' @return matrix of characteristic functions
#'
#' @examples
#' l <- list(c("A","B"), list("C"), list("F", "D"))
#' make.label.table(l)
#'
#' @export make.label.table
make.label.table <- function(labels){
    all_labels <- labels %>% unlist() %>% unique()
    match.mat <- NULL
    for(l in labels){
        match.mat <- rbind(match.mat, find.occ(l,all_labels))
    }
    match.mat
}

#' Find occurrences
#'
#' Prepares a Boolean array representing the characteristic function
#' of a subset. 'some' and 'all' must be sortable. The order of the elements
#' in the output array follow the sorted order of the universe set 'all'
#'
#' @param some subset
#' @param all universe set
#' @return an array representing the characteristic function of a subset of 'all'
#'
#' @examples
#' a <- c("A","B","D")
#' b <- c("A","B","C", "F", "D")
#' find.occ(a,b)
#'
#' @export find.occ
find.occ <- function(some, all){
    some <- sort(some)
    all <- sort(all)
    i <- 1
    j <- 1
    is.match <- rep(FALSE, length(all))
    while( i <= length(some) && j <= length(all)){
        if(some[i] == all[j]){
            is.match[j]<-TRUE
            i <- i+1
        }
        j<-j+1
    }
    names(is.match) <- all
    is.match
}

#' Check if character arrays are equal
#'
#' arrays are sorted beforehand, the procedure is case sensitive
#'
#' @param a first char array
#' @param b second char array
#' @return TRUE if the arrays are equal
#' @examples
#' a <- c("A","B","C")
#' b <- c("A","B","C")
#' c <- c("A","B","D")
#' char.array.equal(a,b)
#' char.array.equal(a,c)
#'
#' @export char.array.equal
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

#' Check if the given label is yet in use
#'
#' Check if the given label is yet in use given a label and a list of labels
#' in which to search
#'
#' @param label a character array of labels or AP
#' @param labels a list of labels
#'
#' @return TRUE if the label is found
#'
#' @examples
#' a <- c("A","B","C")
#' b <- c("A","B","C")
#' c <- c("A","B","D")
#' d <- c(a,b,c)
#' is.not.unique.label(c,d)
#'
#' @export is.not.unique.label
is.not.unique.label <- function(label, labels){
    any(map_lgl(labels ,
                ~ char.array.equal(., label)))
}

#' Find id of node with the given label
#'
#' Check if the given label is yet in use given a label and a list of labels
#' in which to search. Return the first index that matches.
#'
#' @param label a character array of labels or AP
#' @param labels a list of labels
#'
#' @return the first index for which labels matches. 0 if there is none.
#'
#' @examples
#' require(purrr)
#' a <- c("A","B","C")
#' b <- c("A","B","C")
#' c <- c("A","B","D")
#' d <- c(a,b,c)
#' find.label(c,d)
#'
#' @export find.label
find.label <- function(label, labels){
    ans <- which(TRUE == map_lgl(labels ,
                ~ char.array.equal(., label)))
    if(length(ans) == 0){
        0
    }else{
        ans[1]
    }
}

#' Apply function to all cpmc nodes
#'
#' This auxiliary function applies a procedure to each node of the cpmc.
#' Side effects of the procedure are kept during the execution.
#'
#' @param D a cpmc
#' @param procedure a procedure that takes a list in input. This list is
#'        such that the field D contains the cpmc, was.visited is a
#'        Boolean vector containing the nodes that were already
#'        visited by this procedure while current contains the id of the currently
#'        visited node. The output of this procedure must be a named list
#'        containing the updated version of was.visited and current.
#' @export cpmc.apply
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

#' All nodes to sink
#'
#' This function modifies a cpmc so that all nodes point to a sink node
#' (named 'sink_node')
#'
#' this is an example to show the usage of apply
#'
#' @param D a cpmc
#' @export all.to.sink
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

#' Remove singleton
#'
#' This procedure remove singleton (nodes that point only to thyself or nowhere)
#' from a cpmc
#'
#' @param D a cpmc
#' @export cpmc.remove.singleton
cpmc.remove.singleton <- function(D){
    delete.singleton <- function(input){
        i <- input$current
        # remove nodes that do not point to other nodes
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

#' Normalize exiting probabilities
#'
#' This procedure normalizes exiting probabilities for each node of this
#' cpmc. Nothing is done for a node if the sum of these probabilities is 0.
#'
#' @param D a cpmc
#' @export cpmc.normalize
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

#' Apply function to all cpmc nodes that are reachable from a source
#'
#' This auxiliary function applies a procedure to each node that is
#' reachable from a source of the cpmc.
#' Side effects of the procedure are kept during the execution.
#'
#' @param D a cpmc
#' @param procedure a procedure that takes a list in input. This list is
#'        such that the field D contains the cpmc, was.visited is a
#'        Boolean vector containing the nodes that were already
#'        visited by this procedure while current contains the id of the currently
#'        visited node. The output of this procedure must be a named list
#'        containing the updated version of was.visited and current.
#'        (It is actually possible manipulate the DFS visit freely by
#'        changing the output list accordingly)
#' @param id id of the source
#' @param label label of the source (used id is the lowest that matches)
#'
#' @export cpmc.apply.on.visit
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

#' Prune cpmc
#'
#' This procedure removes nodes not reachable form a source
#'
#' @param D a cpmc
#' @param id id of the source
#' @param label label of the source (alternatively)
#'
#' @export cpmc.prune
cpmc.prune <- function(D, id=0, label=NULL){
    simple.visit <- function(input){
        input$was.visited[input$current] <- TRUE
        list(was.visited = input$was.visited)
    }
    out <- cpmc.apply.on.visit(D, simple.visit, id=id, label=label)
    rem <- which(map_lgl(out, ~ !.))
    # very important to remove nodes in descending order!
    walk(sort(rem, decreasing = TRUE), ~ D$remove.node(id = .))
}

#------ TODO --------

#' Apply treatment effects
#'
#' This (experimental) procedure applies treatment effects
#' (sequential)
#'
#' @param D a cpmc
#' @param treatment treatment name
#' @param treatments list of all treatment (as extracted by the lexer-parser)
#' @export cpmc.apply.treatment
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
                        # new edge
                        if(D$P[i,sink] == 0){
                            D$set.edge(i,sink,effect$eff)
                        }else{
                            D$set.edge(i,sink,
                                # p + q - pq
                                effect$eff + D$P[i,sink] - (effect$eff * D$P[i,sink]))
                        }
                    }

                    input$was.visited[i] <- TRUE
                    list(was.visited = input$was.visited, current=i+1)
                }
            }
            apply.treatment <- apply.treatment(effect)

        }else if(effect$action == "mutate"){
            apply.treatment <- function(effect){
                force(effect)
                function(input){
                    i <- input$current
                    if(effect$guard(D$labels[[i]])){
                        # compute new label
                        l <- unlist(setdiff(union(D$labels[[i]], effect$changes$add), effect$changes$remove))

                        # add new node
                        if(D$get.id(l) == 0){
                            # manage "clonal" case TODO make it more general
                            if(is.null(l)){
                                l <- c("Clonal")
                            }
                            D$add.node(l)
                            dest <- D$get.id(l)
                            if(sum(D$P[dest,]) == 0){
                                D$set.edge(dest, dest, 1, normalize = F)
                            }
                            input$was.visited <- c(input$was.visited, FALSE)
                        }
                        dest <- D$get.id(l)

                        # add arc to dest
                        if(D$P[i,dest] == 0){
                            D$set.edge(i,dest,effect$eff)
                        }else{
                            D$set.edge(i,dest,
                                       # p + q - pq
                                       effect$eff + D$P[i,dest] - (effect$eff * D$P[i,dest]))
                        }

                    }

                    input$was.visited[i] <- TRUE
                    list(was.visited = input$was.visited, current=i+1)
                }
            }
            apply.treatment <- apply.treatment(effect)

        }else{
            print(paste("ERROR: mode not implemented:", effect$action))
            return(FALSE)
        }

        cpmc.apply(D, apply.treatment)
    }

    return(TRUE)
}




