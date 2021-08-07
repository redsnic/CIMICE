


#' Create a stub for a generator 
#'
#' Create a generator topology directly from a dataset.
#' The topology will follow the subset relation.
#' 
#' @param dataset A compacted CIMICE dataset
#'
#' @return a generator, with weight = 0 for all the edges
#' 
#' @examples 
#' make_generator_stub(example_dataset())
#' 
#' @export make_generator_stub
make_generator_stub <- function(dataset){
  graph_data <- example_dataset() %>% dataset_preprocessing()
  generator <- graph_non_transitive_subset_topology(graph_data$samples, graph_data$labels) %>% 
    as_tbl_graph() %>% 
    mutate(matching_sample = graph_data$matching_samples)
  for(label in generator %>% pull(label) ){
    generator <- generator %>% bind_edges(from = label, to = label, node_key = "label")
  }
  
  generator <- generator %E>% mutate(weight = 0)
  generator <- activate(generator, what="nodes")
  generator
}

#' Add edge weights to a generator
#' 
#' @param generator a generator 
#' @param f_t_w_list a list of triplets (from, to, list), the triplets must not 
#' be nested in the list. For example list("A","B",0.3, "B", "C", 0.2) is a valid
#' input.
#' @param by "labels" or "samples" to use gene labels or sample labels as references for
#' edge identifiers.
#' 
#' 
#' @return the generator with the modified edges (invalid edges are ignored)
#'   
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'  make_generator_stub() %>% 
#'  set_generator_edges(
#'  list(
#'    "S1", "S4", 1 , 
#'    "S2, S3", "S4", 1 , 
#'    "S4", "S7", 1 , 
#'    "S4", "S5, S6, S8", 1 , 
#'    "Clonal", "S1", 1 , 
#'    "Clonal", "S2, S3", 1 , 
#'    "S1", "S1", 1 , 
#'    "S2, S3", "S2, S3", 1 , 
#'    "S4", "S4", 1 , 
#'    "S7", "S7", 1 , 
#'    "S5, S6, S8", "S5, S6, S8", 1 , 
#'    "Clonal", "Clonal", 1 
#'  ))
#'    
#' @export set_generator_edges      
set_generator_edges <- function(generator, f_t_w_list, by = "labels"){
  
  # create columns for the data frame
  type <- 1
  from_l <- list()
  to_l <- list()
  weight_l <- c()
  for(el in f_t_w_list){
    if(type == 1){
      from_l <- c(from_l, el)
      type <- 2
    }else if(type == 2){
      to_l <- c(to_l, el)
      type <- 3
    }else{
      assert_that(is.numeric(el), msg="weights must be numeric")
      weight_l <- c(weight_l, el)
      type <- 1
    }
  }
  assert_that(length(from_l) == length(to_l) & length(to_l) == length(weight_l), 
              msg="The sequence must be multiple of three (from, to, weight)")
  
  # get node ids by label/sample_label
  if(by == "samples"){
    id_col <- "matching_samples"
  }else{
    id_col <- "label"
  }
  
  from_l <- map_int(from_l, ~ ifelse(is.numeric(.),
                                     ., 
                                     which((generator %N>% pull(!!id_col)) == .)))
  to_l <- map_int(to_l, ~ ifelse(is.numeric(.),
                                     ., 
                                     which((generator %N>% pull(!!id_col)) == .)))
  
  # create the weight dataframe
  weights_df <- data.frame(from=from_l,to=to_l,weight=weight_l)
  print(weights_df)
  
  # add weights
  generator <- generator %E>% left_join(weights_df, by=c("from","to")) %>% 
    mutate(weight = ifelse(is.na(weight.y), weight.x, weight.y)) %>% 
    select(-weight.x, -weight.y)
  
  generator <- activate(generator, what="nodes")
  generator
}

#' Prepare a command to add edge weights to a generator
#' 
#' Prints a string in the form of the command that 
#' sets weights for all the edges of this generator.
#'
#' @param generator a generator 
#' @param by "labels" or "samples" to use gene labels or sample labels as references for
#' edge identifiers.
#'
#' @examples 
#' require(dplyr)
#' example_dataset() %>% 
#'   make_generator_stub() %>%
#'   prepare_generator_edge_set_command()
#'
#' @export prepare_generator_edge_set_command
prepare_generator_edge_set_command <- function(generator, by="labels"){
  edges <- generator %E>% data.frame()
  generator %N>% print
  if(by == "samples"){
    names_vec <- generator %N>% pull(matching_sample)
  }else{
    names_vec <- generator %N>% pull(label)
  }
  out <- c("set_generator_edges(\n  list(\n")
  for(row in seq(1,nrow(edges))){
    out <- c(out, "    ")
    out <- c(out, paste(
      paste('"', names_vec[edges[row, "from"]], '"', sep=""),
      paste('"', names_vec[edges[row, "to"]], '"', sep=""),
      "_", sep=", "))
    if(row < nrow(edges)){
      out <- c(out, ",")
    }
    out <- c(out, "\n")
  }
  out <- c(out, "))\n")
  paste(out) %>% cat
}

#' Plot a generator
#'
#' Simple ggraph interface to draw a generator 
#' 
#' @param generator a generator
#'
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'  make_generator_stub() %>% 
#'  set_generator_edges(
#'  list(
#'    "S1", "S4", 1 , 
#'    "S2, S3", "S4", 1 , 
#'    "S4", "S7", 1 , 
#'    "S4", "S5, S6, S8", 1 , 
#'    "Clonal", "S1", 1 , 
#'    "Clonal", "S2, S3", 1 , 
#'    "S1", "S1", 1 , 
#'    "S2, S3", "S2, S3", 1 , 
#'    "S4", "S4", 1 , 
#'    "S7", "S7", 1 , 
#'    "S5, S6, S8", "S5, S6, S8", 1 , 
#'    "Clonal", "Clonal", 1 
#'  )) %>% finalize_generator %>%  plot_generator
#'  
#' @export plot_generator
plot_generator <- function(generator){
  ggraph(generator, layout = "sugiyama") +
    geom_node_point() +
    geom_node_label(aes(label=label)) +
    geom_edge_link(
      aes(label=weight, start_cap = label_rect(node1.label),
          end_cap = label_rect(node2.label)),
      arrow = arrow(length = unit(4, 'mm')),
      label_dodge = unit(2.5, 'mm'),
      angle_calc = 'along'
    ) +
    geom_edge_loop(
      aes(label=weight, start_cap = label_rect(node1.label),
          end_cap = label_rect(node2.label)),
      arrow = arrow(length = unit(4, 'mm')),
      label_dodge = unit(2.5, 'mm'),
      angle_calc = 'along'
    ) 
}

#' Finalize generator normalizing edge weights
#' 
#' Checks if a generator can be normalized so that it 
#' actually is a Markov Chain 
#' 
#' @param generator a generator 
#' 
#' @examples 
#' require(dplyr)
#' 
#' example_dataset() %>%
#'  make_generator_stub() %>% 
#'  set_generator_edges(
#'  list(
#'    "S1", "S4", 1 , 
#'    "S2, S3", "S4", 1 , 
#'    "S4", "S7", 1 , 
#'    "S4", "S5, S6, S8", 1 , 
#'    "Clonal", "S1", 1 , 
#'    "Clonal", "S2, S3", 1 , 
#'    "S1", "S1", 1 , 
#'    "S2, S3", "S2, S3", 1 , 
#'    "S4", "S4", 1 , 
#'    "S7", "S7", 1 , 
#'    "S5, S6, S8", "S5, S6, S8", 1 , 
#'    "Clonal", "Clonal", 1 
#'  )) %>% finalize_generator
#' 
#' @export finalize_generator
finalize_generator <- function(generator){
  normalization_df <- generator %E>% 
    data.frame() %>% 
    group_by(from) %>% 
    summarise(tot = sum(weight)) 
  
  print(normalization_df)
  
  if(length(which(normalization_df %>% pull(tot) == 0)) != 0){
    stop("Some node has no exiting edge with probability greater than 0.\nCheck your call to set_generator_edges")
  }
  
  generator %E>% 
    left_join(normalization_df, by="from") %>% 
    mutate(weight = weight/tot) %>% 
    select(-tot)
}



