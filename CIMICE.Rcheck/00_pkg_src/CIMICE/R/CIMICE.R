#' CIMICE Package
#'
#' CIMICE-R: (Markov) Chain Method to Inferr Cancer Evolution
#'
#' @description R implementation of the CIMICE tool
#'
#' @author Nicolò Rossi \email{olocin.issor@gmail.com}
#'
#' @docType package
#' @name CIMICE
#'
#' @importFrom dplyr select enexprs select_if %>% ungroup
#' @importFrom ggplot2 ggplot geom_histogram aes arrow unit labs
#' @importFrom glue glue
#' @importFrom tidyr drop_na
#' @importFrom igraph V E as_ids tail_of head_of get.edge.ids
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
#' @importFrom igraph graph_from_edgelist as_adj E<- V<-
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom visNetwork toVisNetworkData visNetwork
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom purrr map2 map_dbl
#' @importFrom ggraph ggraph geom_node_point geom_node_label
#' @importFrom ggraph geom_edge_link label_rect
#' @importFrom plyr count
#' @importFrom stats cor
#' @importFrom utils read.csv
#' @importFrom relations transitive_reduction endorelation
#'
NULL
