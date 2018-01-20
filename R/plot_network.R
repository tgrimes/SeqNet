# Requires: 
#   igraph - for plotting networks.
# require(igraph)

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
#' @return Creates a plot of the network and returns a graph object. 
#' The graph object can be passed back into a future call of plot.network() 
#' through the `compare_edge` argument, which will setup the plot for easier 
#' comparison between the old graph and the graph of `network`.
#' @export
plot.network <- function(network, compare_graph = NULL,
                         weighted = NULL, as_subgraph = FALSE,
                         node_scale = 5, edge_scale = 1, 
                         coords =  igraph::layout.fruchterman.reingold,
                         main = "Untitled", ...) {
  library(igraph)
  
  if(class(network) == "network") {
    adj_matrix <- get_adj_matrix_from_network(network)
  } else if(is.matrix(network)){
    if(!all(network %in% c(0, 1))) {
      network[network != 0] <- 1 # Set any nonzero values to 1.
    }
    adj_matrix <- network
  } else {
    stop("network should be either a matrix or a network class object.")
  }
  
  if(is.null(colnames(adj_matrix))) {
    colnames(adj_matrix) <- 1:ncol(adj_matrix)
  }
  
  g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                           mode = "undirected",
                                           weighted = weighted)
  igraph::V(g)$size <- log(apply(adj_matrix, 2, sum) + 1) 
  igraph::V(g)$size <- igraph::V(g)$size / max(igraph::V(g)$size) * node_scale
  igraph::V(g)$frame.color <- "white"
  
  if(sum(adj_matrix) == 0 & !is.null(compare_graph)) {
    # The given network has no edges; plot the reference network with no edges.
    vertex.label.color <- rgb(0, 0, 0, 0) # Make labels transparent.
    
    plot(compare_graph, vertex.color = "white", vertex.label.font = 2,
         vertex.label.color = vertex.label.color, vertex.label.cex = 0.7,
         edge.color = "wheat", main = main, ...)
  } else {
    igraph::E(g)$width <- edge_scale
    if(as_subgraph) {
      if(!is.null(compare_graph)) {
        g <- igraph::induced_subgraph(g, union(which(igraph::V(g)$size > 0), 
                                      which(igraph::V(compare_graph)$size > 0)))
      } else {
        g <- igraph::induced_subgraph(g, which(igraph::V(g)$size > 0))
      }
    }
    
    edge.color <- "black"
    vertex.label.color = "blue"
    if(!is.null(compare_graph)) {
      h <- compare_graph %u% g
      edge.color <- rep("red", length(igraph::E(h))) # Default color for edge in compared network is red.
      subset_1 <- which(attr(igraph::E(h), "vnames") %in% attr(igraph::E(compare_graph), "vnames"))
      edge.color[subset_1] <- 
        ifelse(attr(igraph::E(h), "vnames")[subset_1] %in% attr(igraph::E(g), "vnames"), 
               rep("black", length(igraph::E(compare_graph))),
               rep("wheat", length(igraph::E(compare_graph))))
      g <- h
      
      vertex.label.color <- rgb(0, 0, 0, 0) # Make labels transparent.
      
      igraph::V(g)$size <- log(apply(adj_matrix, 2, sum) + 1) 
      igraph::V(g)$size <- igraph::V(g)$size / max(igraph::V(g)$size) * node_scale
      igraph::V(g)$frame.color <- "white"
      
      coords <- igraph::layout.fruchterman.reingold(compare_graph)
    }
    
    plot(g, vertex.color = "orange", vertex.label.font = 2,
         vertex.label.color = vertex.label.color, vertex.label.cex = 0.7,
         edge.color = edge.color, layout = coords,
         main = main, ...)
  }
  
  return(g)
}

#' Plot function for networks
#' 
#' This function creates a plot for a network object or adjacency matrix. 
#' If the result of another plot is provided, this plot will be modified for 
#' easier comparison.
#' @param network Either a network object or adjacency matrix of the network.
#' @param compare_graph The plot of another network to use for comparison.
#' @param weighted Are the edges weighted? 
#' @param as_subgraph If true, only nodes of positive degree will be shown.
#' @param node_scale Used for scaling of nodes.
#' @param edge_scale used for scaling of edges.
#' @param coords Layout used for the network.
#' @return Creates a plot of the network and returns a graph object. 
#' The graph object can be passed back into a future call of plot.network() 
#' through the `compare_edge` argument, which will setup the plot for easier 
#' comparison between the old graph and the graph of `network`.
#' @export
plot_network <- plot.network
