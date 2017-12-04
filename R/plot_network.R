# Requires: 
#   igraph - for plotting networks.
library(igraph)

#' Plot function for network class
#' 
#' This function plots the given network. If the result of another plot is 
#' provided, this plot will be modified for easier comparison.
#' @param network Either a network object or adjacency matrix of the network.
#' @param compare_graph The plot of another network to use for comparison.
#' @param weighted Are the edges weighted? 
#' @param as_subgraph If true, only nodes of positive degree will be shown.
#' @param node_scale Used for scaling of nodes.
#' @param edge_scale used for scaling of edges.
#' @param coords Layout used for the network.
#' @value Creates a plot of the network and returns a graph object. 
#' The graph object can be passed back into a future call of plot.network() 
#' through the `compare_edge` argument, which will setup the plot for easier 
#' comparison between the old graph and the graph of `network`.
#' @export
plot.network <- function(network, compare_graph = NULL,
                         weighted = NULL, as_subgraph = FALSE,
                         node_scale = 5, edge_scale = 1, 
                         coords =  layout.fruchterman.reingold,
                         main = "Untitled", ...) {
  if(class(network) == "network") {
    adj_matrix <- get_adj_matrix_from_network(network)
  } else {
    adj_matrix <- network
  }
  
  if(is.null(colnames(adj_matrix))) {
    colnames(adj_matrix) <- 1:ncol(adj_matrix)
  }
  
  g <- graph_from_adjacency_matrix(adj_matrix,
                                   mode = "undirected",
                                   weighted = weighted)
  V(g)$size <- log(apply(adj_matrix, 2, sum) + 1) 
  V(g)$size <- V(g)$size / max(V(g)$size) * node_scale
  V(g)$frame.color <- "white"
  
  if(sum(adj_matrix) == 0 & !is.null(compare_graph)) {
    plot(compare_graph, vertex.color = "white", vertex.label.font = 2,
         vertex.label.color = "blue", vertex.label.cex = 0.7,
         edge.color = "white", main = main, ...)
  } else {
    E(g)$width <- edge_scale
    if(as_subgraph) {
      if(!is.null(compare_graph)) {
        g <- induced_subgraph(g, union(which(V(g)$size > 0), 
                                       which(V(compare_graph)$size > 0)))
      } else {
        g <- induced_subgraph(g, which(V(g)$size > 0))
      }
    }
    
    edge.color <- "black"
    if(!is.null(compare_graph)) {
      h <- compare_graph %u% g
      edge.color <- rep("red", length(E(h))) # Default color for edge in compared network is red.
      subset_1 <- which(attr(E(h), "vnames") %in% attr(E(compare_graph), "vnames"))
      edge.color[subset_1] <- 
        ifelse(attr(E(h), "vnames")[subset_1] %in% attr(E(g), "vnames"), 
               rep("black", length(E(compare_graph))),
               rep("wheat", length(E(compare_graph))))
      g <- h
      
      V(g)$size <- log(apply(adj_matrix, 2, sum) + 1) 
      V(g)$size <- V(g)$size / max(V(g)$size) * node_scale
      V(g)$frame.color <- "white"
      
      if(!is.null(compare_graph)) {
        coords <- layout.fruchterman.reingold(compare_graph)
      }
    }
    
    plot(g, vertex.color = "orange", vertex.label.font = 2,
         vertex.label.color = "blue", vertex.label.cex = 0.7,
         edge.color = edge.color, layout = coords,
         main = main, ...)
  }
  
  return(g)
}
