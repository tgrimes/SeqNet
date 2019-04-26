# Requires: 
#   igraph - for plotting networks.
# require(igraph)

setClass(Class = "network_plot")

#' Visualize a network
#' 
#' This function is used to plot a network. The 'network' argument can be a 
#' network object, network module, an adjacency matrix, or an association matrix. 
#' If the result of another plot is provided using the 'compare_graph' argument, 
#' then the layout of this network will be based on that plot.
#' @param network A 'network', 'network_module', or 'matrix' object.
#' @param compare_graph The plot of another network to use for comparison.
#' @param as_subgraph If TRUE, only nodes of positive degree will be shown.
#' @param node_scale Used for scaling of nodes.
#' @param edge_scale Used for scaling of edges.
#' @param node_color The color used for the nodes.
#' @param generate_coords A function to generate the layout of a graph; used
#' if coords is NULL. See ?igraph::layout_ for details. Other options include 
#' 'igraph::as_star', 'igraph::in_circle', and 'igraph::with_fr', among many others.
#' @param include_vertex_labels If TRUE, the verticies will be labeled.
#' @param display_plot If TRUE (default), the plot will be generated and displayed.
#' @param ... Additional arguments passed to plot.igraph().
#' @return Creates a plot of the network and returns a graph object. 
#' The graph object can be passed back into a future call of plot.network() 
#' through the `compare_edge` argument, which will setup the plot for easier 
#' comparison between the old graph and the graph of `network`.
#' @export
plot_network <- function(network, compare_graph = NULL, as_subgraph = FALSE,
                         node_scale = 4, edge_scale = 1,
                         node_color = adjustcolor("orange", 0.5),
                         generate_layout = igraph::nicely,
                         include_vertex_labels = TRUE, 
                         display_plot = TRUE, ...) {
  ##################################
  # Check arguments for errors.
  
  if(!(class(network) %in% c("network", "network_module", "matrix"))) {
    stop(paste0("Argument 'network' must be a 'network', 'network_module', ",
                "or 'matrix' object."))
  }
  
  if((class(network) == "matrix") && is.null(colnames(network))) {
    # Set default gene names in a matrix to 1:p.
    colnames(network) <- 1:ncol(network)
  }

  if(!is.null(compare_graph)) {
    if(class(compare_graph) != "network_plot")
      stop("Argument 'compare_graph' must be a 'network_plot' object.")
    if((nrow(compare_graph$coords) != length(get_node_names(network))) &&
       !all(attr(igraph::V(compare_graph$graph), "names") %in% get_node_names(network)))
      stop(paste("Argument 'compare_graph' and 'network' must contain the same",
                 "number of nodes or contain an overlapping set of nodes."))
  }
  
  if(node_scale <= 0) 
    stop("Argument 'node_scale' must be positive.")
  if(edge_scale <= 0) 
    stop("Argument 'edge_scale' must be positive.")
  
  ##################################
  
  # Initialize plot and obtain an association matrix if the network is weighted.
  if(is_weighted(network)) {
    assoc_matrix <- abs(get_association_matrix(network))
    assoc_matrix[abs(assoc_matrix) < 10^-13] <- 0 # Set small associations to zero.
    adj_matrix <- 1 * (assoc_matrix != 0)
    g <- igraph::graph_from_adjacency_matrix(assoc_matrix,
                                             mode = "undirected",
                                             weighted = TRUE)
  } else {
    adj_matrix <- get_adjacency_matrix(network)
    assoc_matrix <- NULL
    g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                             mode = "undirected",
                                             weighted = NULL)
  }
  
  # Initialize the comparison plot, if one is provided.
  if(is.null(compare_graph)) {
    g_compare <- NULL
  } else {
    # Initialize the graph to compare to. 
    g_compare <- compare_graph$graph
    # Subset both graphs to common nodes.
    common_nodes <- intersect(attr(igraph::V(g), "names"), attr(igraph::V(g_compare), "names"))
    index_subset_g <- which(attr(igraph::V(g), "names") %in% common_nodes)
    index_subset_g_compare <- which(attr(igraph::V(g_compare), "names") %in% common_nodes)
    
    # Update based on comparison graph.
    g <- igraph::induced_subgraph(g, index_subset_g)
    if(!is.null(assoc_matrix)) {
      assoc_matrix <- assoc_matrix[index_subset_g, index_subset_g]
    }
    adj_matrix <- adj_matrix[index_subset_g, index_subset_g]
  }
  
  # Create subgraph, if requested.
  if(as_subgraph) {
    # Determine which nodes to keep in subgraph.
    index_nodes <- unname(which(igraph::degree(g) > 0))
    
    # If g contains no edges, then a subgraph cannot be created.
    if(length(index_nodes) == 0) {
      warning("Cannot create subgraph; the network must contain at least one edge.")
    } else {
      # Update based on subset of nodes.
      if(!is.null(g_compare)) {
        index_subset_g_compare <- index_subset_g_compare[index_nodes]
      }
      g <- igraph::induced_subgraph(g, index_nodes)
      if(!is.null(assoc_matrix)) {
        assoc_matrix <- assoc_matrix[index_nodes, index_nodes]
      }
      adj_matrix <- adj_matrix[index_nodes, index_nodes]
    }
  }
  
  # Initialize coordinates for graph layout.
  if(is.null(compare_graph)) {
    coords = igraph::layout_(g, generate_layout())
  } else {
    coords = compare_graph$coords[index_subset_g_compare, ]
  }
  if(nrow(coords) != igraph::vcount(g)) {
    stop("coords do not match number of verticies in the graph.")
  }
  
  # Adjust node size, color, and frame.
  if(!is.null(assoc_matrix)) {
    # Scale associations relative to largest association in the network.
    temp <- max(assoc_matrix[lower.tri(assoc_matrix)])
    node_weights <- sqrt(apply(assoc_matrix, 2, sum) / 
                           ifelse(temp == 0, 1, temp))
  } else {
    # Default node weights are proportional to sqrt(degree)
    node_weights <- sqrt(igraph::degree(g))
  }
  vertex.size <- node_weights * node_scale
  vertex.frame.color <- "white"
  if(include_vertex_labels) {
    vertex.label.color <- rgb(0.1, 0.1, 0.1, 0.8)
  } else {
    vertex.label.color <- rgb(0, 0, 0, 0) # Make labels transparent.
  }
  
  # Adjust edge width and color.
  if(!is.null(assoc_matrix)) {
    edge_weights <- assoc_matrix[lower.tri(assoc_matrix)]
    edge_weights <- edge_weights[edge_weights != 0]
    edge_weights <- edge_weights / 
      ifelse(max(edge_weights) == 0, 1, max(edge_weights))
    if(length(edge_weights) != length(igraph::E(g))) {
      stop("Edge weights do not match number of edges.")
    }
  } else {
    edge_weights <- rep(1, length(igraph::E(g)))
  } 
  edge.width <- edge_weights * edge_scale
  edge.color <- "black"
  
  if(display_plot) {
    plot(g, vertex.color = node_color, vertex.label.font = 2,
         vertex.size = vertex.size, vertex.frame.color = vertex.frame.color,
         vertex.label.color = vertex.label.color, vertex.label.cex = 0.7,
         edge.color = edge.color, edge.width = edge.width, layout = coords, 
         ...)
  }
  
  plot_summary <- list(graph = g,
                       coords = coords,
                       modules = NULL,
                       node_weights = node_weights,
                       edge_weights = edge_weights,
                       vertex.size = vertex.size,
                       vertex.frame.color = vertex.frame.color,
                       vertex.label.color = vertex.label.color,
                       vertex.color = node_color,
                       edge.color = edge.color,
                       edge.width = edge.width)
  class(plot_summary) <- "network_plot"
  
  # Return list silently:
  invisible(plot_summary)
}



#' Visualize a network and its modules
#' 
#' This function plots a network and highlights the individual modules. 
#' An attempt is made to layout the nodes so that any visual overlaps among modules
#' correspond to true overlaps in the network, however it is possible that 
#' a node may appear to be in multiple modules in the visualization when it does
#' not actually belong to multiple modules. If the result of another plot is 
#' provided using the 'compare_graph' argument, then the layout of this network
#' will be based on that plot and convex hulls are drawn to trace out the modules; 
#' in this case it is likely that the displayed modules will contain extraneous
#' nodes.
#' @param network A 'network' object to plot. Alternatively, an adjacency or
#' association matrix can be provided, in which case the 'modules' argument
#' should be specified.
#' @param compare_graph The plot of another network to use for comparison.
#' @param as_subgraph If TRUE, only nodes of positive degree will be shown.
#' @param modules A list of modules for the network; this is used to provide
#' a member list of each module when the 'network' argument is not a 'network' 
#' object. To get this list from a network, use 'get_network_modules()'.
#' @param node_scale Used for scaling of nodes.
#' @param edge_scale Used for scaling of edges.
#' @param node_color The color used for the nodes.
#' @param group_color A vector of colors used for the modules.
#' @param generate_coords A function to generate the layout of a graph; used
#' if coords is NULL. See ?igraph::layout_ for details. Other options include 
#' 'igraph::as_star', 'igraph::in_circle', and 'igraph::with_fr', among many others.
#' @param include_vertex_labels If TRUE, the verticies will be labeled.
#' @param show_legend If TRUE, a legend for the modules is shown. Default is FALSE.
#' @param legend_position The location of the legend. Can be any one of "bottomright",
#' "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" or "center".
#' @param legend_horizontal If TRUE, the legend will be displayed horizontally.
#' @param display_plot If TRUE (default), the plot will be generated and displayed.
#' @param ... Additional arguments passed to plot.igraph().
#' @return A 'network_plot' object for the network. This object can be passed 
#' back into a future call of plot.network() through the `compare_graph` 
#' argument, which will setup the plot for easier comparison between the old 
#' graph and the new graph of `network`.
#' @export
plot_modules <- function(network, compare_graph = NULL, as_subgraph = TRUE,
                         modules = NULL,
                         node_scale = 4, edge_scale = 1,
                         node_color = adjustcolor("orange", 0.5),
                         group_color = RColorBrewer::brewer.pal(9, 'Set1'),
                         coords =  NULL, 
                         generate_layout = igraph::nicely,
                         include_vertex_labels = TRUE, 
                         show_legend = FALSE,
                         legend_position = "topright", 
                         legend_horizontal = FALSE, 
                         display_plot = TRUE, ...) {
  ##################################
  # Check arguments for errors.
  
  if(!(class(network) %in% c("network", "matrix"))) {
    stop("Argument 'network' must be a 'network' object.")
  }
  if((class(network) %in% "matrix") && is.null(modules)) {
    stop("If 'network' is a matrix, then the 'modules' argument must be specified.")
  }
  
  if(is.null(modules)) {
    modules <- lapply(network$modules, function(module) module$nodes)
  }
  
  # Module fill gets an additional alpha value for transparency:
  group_color_fill <- paste0(group_color, '10')
  
  node_names <- get_node_names(network)
  p <- length(node_names)
  
  if(!is.null(compare_graph)) {
    if(class(compare_graph) != "network_plot")
      stop("Argument 'compare_graph' must be a 'network_plot' object.")
    if((nrow(compare_graph$coords) != p) &&
       !all(attr(igraph::V(compare_graph$graph), "names") %in% node_names))
      stop(paste("Argument 'compare_graph' and 'network' must contain the same",
                 "number of nodes or contain an overlapping set of nodes."))
  }
  
  if(node_scale <= 0) 
    stop("'node_scale' must be positive.")
  if(edge_scale <= 0) 
    stop("'edge_scale' must be positive.")
  
  ##################################
  
  # Initialize plot and obtain an association matrix if the network is weighted.
  nodes <- 1:p # Used for updating modules.
  if(is_weighted(network)) {
    assoc_matrix <- abs(get_association_matrix(network))
    assoc_matrix[abs(assoc_matrix) < 10^-13] <- 0 # Set small associations to zero.
    adj_matrix <- 1 * (assoc_matrix != 0)
    g <- igraph::graph_from_adjacency_matrix(assoc_matrix,
                                             mode = "undirected",
                                             weighted = TRUE)
  } else {
    adj_matrix <- get_adjacency_matrix(network)
    assoc_matrix <- NULL
    g <- igraph::graph_from_adjacency_matrix(adj_matrix,
                                             mode = "undirected",
                                             weighted = NULL)
  }

  # Initialize the comparison plot, if one is provided.
  if(is.null(compare_graph)) {
    g_compare <- NULL
  } else {
    # Initialize the graph to compare to. 
    g_compare <- compare_graph$graph
    # Subset both graphs to common nodes.
    common_nodes <- intersect(attr(igraph::V(g), "names"), attr(igraph::V(g_compare), "names"))
    index_subset_g <- which(attr(igraph::V(g), "names") %in% common_nodes)
    index_subset_g_compare <- which(attr(igraph::V(g_compare), "names") %in% common_nodes)
    
    # Update based on comparison graph.
    g <- igraph::induced_subgraph(g, index_subset_g)
    nodes <- nodes[index_subset_g]
    modules <- lapply(modules, function(m) m[m %in% nodes])
    # If any modules become empty, remove them from the list.
    modules <- modules[sapply(modules, function(m) length(m) > 0)]
    if(!is.null(assoc_matrix)) {
      assoc_matrix <- assoc_matrix[index_subset_g, index_subset_g]
    }
    adj_matrix <- adj_matrix[index_subset_g, index_subset_g]
  }
  
  # Create subgraph, if requested.
  if(as_subgraph) {
    # Determine which nodes to keep in subgraph.
    index_nodes <- unname(which(igraph::degree(g) > 0))

    # If g contains no edges, then a subgraph cannot be created.
    if(length(index_nodes) == 0) {
      warning("Cannot create subgraph; the network must contain at least one edge.")
    } else {
      # Update based on subset of nodes.
      if(!is.null(g_compare)) {
        index_subset_g_compare <- index_subset_g_compare[index_nodes]
      }
      g <- igraph::induced_subgraph(g, index_nodes)
      nodes <- nodes[index_nodes]
      modules <- lapply(modules, function(m) m[m %in% nodes])
      # If a module had no connections, all of its modules will be removed. In
      # this case, remove that module from the list.
      modules <- modules[sapply(modules, function(m) length(m) > 0)]
      if(!is.null(assoc_matrix)) {
        assoc_matrix <- assoc_matrix[index_nodes, index_nodes]
      }
      adj_matrix <- adj_matrix[index_nodes, index_nodes]
    }
    
  }
  
  # Update modules to use subgraph node labels.
  modules <- lapply(modules, function(m) {
    unlist(sapply(m, function(node) which(nodes == node)))
  })
  
  # Initialize coordinates for graph layout.
  if(is.null(coords)) {
    if(is.null(compare_graph)) {
      coords = get_layout_for_modules(g, modules)
    } else {
      coords = compare_graph$coords[index_subset_g_compare, ]
    }
  }
  if(nrow(coords) != igraph::vcount(g)) {
    stop("coords do not match number of verticies in the graph.")
  }
  
  # Adjust node size, color, and frame.
  if(!is.null(assoc_matrix)) {
    # Scale associations relative to largest association in the network.
    temp <- max(assoc_matrix[lower.tri(assoc_matrix)])
    node_weights <- sqrt(apply(assoc_matrix, 2, sum) / 
                           ifelse(temp == 0, 1, temp))
  } else {
    # Default node weights are proportional to sqrt(degree)
    node_weights <- sqrt(igraph::degree(g))
  }
  vertex.size <- node_weights * node_scale
  vertex.frame.color <- "white"
  if(include_vertex_labels) {
    vertex.label.color <- rgb(0.1, 0.1, 0.1, 0.8)
  } else {
    vertex.label.color <- rgb(0, 0, 0, 0) # Make labels transparent.
  }
  
  # Adjust edge width and color.
  if(!is.null(assoc_matrix)) {
    edge_weights <- assoc_matrix[lower.tri(assoc_matrix)]
    edge_weights <- edge_weights[edge_weights != 0]
    edge_weights <- edge_weights / 
      ifelse(max(edge_weights) == 0, 1, max(edge_weights))
    if(length(edge_weights) != length(igraph::E(g))) {
      stop("Edge weights do not match number of edges.")
    }
  } else {
    edge_weights <- rep(1, length(igraph::E(g)))
  } 
  edge.width <- edge_weights * edge_scale
  edge.color <- "black"
  
  module_names <- names(modules)
  if(is.null(module_names)) {
    module_names <- paste("module", 1:length(modules), sep = "_")
  }
  
  if(display_plot) {
    plot(g, vertex.color = node_color, vertex.label.font = 2,
         vertex.size = vertex.size, vertex.frame.color = vertex.frame.color,
         vertex.label.color = vertex.label.color, vertex.label.cex = 0.75,
         edge.color = edge.color, edge.width = edge.width, layout = coords, 
         mark.groups = modules,
         mark.shape = 1,
         mark.col = group_color_fill,
         mark.border = group_color, ...)
    if(show_legend) {
      legend(legend_position,
             legend = module_names, 
             col = paste0(group_color, '75'),
             pch = 15, bty = "n",  pt.cex = 1.5, cex = 0.7, 
             text.col = "black", horiz = legend_horizontal)
    }
  }
  
  plot_summary <- list(graph = g,
                       coords = coords,
                       modules = modules,
                       node_weights = node_weights,
                       edge_weights = edge_weights,
                       vertex.size = vertex.size,
                       vertex.frame.color = vertex.frame.color,
                       vertex.label.color = vertex.label.color,
                       vertex.color = node_color,
                       edge.color = edge.color,
                       edge.width = edge.width)
  class(plot_summary) <- "network_plot"
  
  # Return list silently:
  invisible(plot_summary)
}

#' Creates coordinates based on a set of modules
#' 
#' @param g An 'igraph' object
#' @param modules A list containing sets of indicies indicating the nodes g that 
#' belong to each module
#' @return A matrix of coordinates for plotting
#' @export
get_layout_for_modules <- function(g, modules) {
  # See the R blog post by Markus Konrad for details:
  # https://www.r-bloggers.com/visualizing-graphs-with-overlapping-node-groups/
  
  # Get node set from graph.
  nodes <- 1:length(igraph::V(g))
  node_labels <- attr(igraph::V(g), "names")
  # Remove any empty modules.
  modules <- modules[sapply(modules, length) > 0]
  
  n_groups <- length(modules)
  if(is.null(names(modules)))
    names(modules) <- 1:n_groups
  modules_df <- data.frame(
    id = unlist(modules),
    group = rep(1:n_groups, sapply(modules, length))
  )
  
  edges_network <- do.call(rbind, 
                           lapply(lapply(attr(igraph::E(g), "vnames"), 
                                         strsplit, split = "\\|"), 
                                  function(x) c(which(node_labels == x[[1]][1]),
                                                which(node_labels == x[[1]][2]))))
  
  memberless_nodes <- setdiff(nodes, unique(unlist(modules)))
  if(length(memberless_nodes) > 0) {
    modules_df <- rbind(modules_df, 
                        data.frame(id = memberless_nodes, group = NA))
  }
  
  group_nodes <- max(nodes) + 1:n_groups
  edges <- 
    purrr::map_dfr(split(modules_df, modules_df$group), function (grp) {
      group_id <- grp$group[1]
      edges_module <- edges_network[(edges_network[, 1] %in% grp$id) &
                                      (edges_network[, 2] %in% grp$id), ]
      # Connect nodes within module.
      df <- data.frame(from = edges_module[, 1], to = edges_module[, 2],
                       weight = 1, group = group_id)
      # Tie the module together.
      from_group <- rep(group_nodes[group_id], nrow(grp))
      to_nodes <- grp$id
      if(length(to_nodes) > 0) {
        df <- rbind(df,
                    data.frame(from = from_group, to = to_nodes,
                               weight = 5,  group = group_id))
      }
      
      return(df)
    })
  # Repel nodes that are in different modules.
  for(node in nodes) {
    outside_nodes <- unique(unlist(modules[!sapply(modules,
                                                   function(module_nodes) node %in% module_nodes)]))
    if(length(outside_nodes) > 0) {
      temp <- rbind(expand.grid(node, outside_nodes),
                    expand.grid(outside_nodes, node))
      temp <- temp[temp[, 1] < temp[, 2], ]
      edges <- rbind(edges,
                     data.frame(from = temp[, 1], to = temp[, 2],
                                weight = 0.001, group = NA))
    }
  }
  edges <- edges[order(edges[, 1]), ]
  nodes_df <- data.frame(id = c(nodes, group_nodes))
  g_virt <- igraph::graph_from_data_frame(edges, 
                                          directed = FALSE,
                                          vertices = nodes_df)
  coords <- igraph::layout_nicely(g_virt)
  coords <- coords[1:length(nodes), ]
  return(coords)
}



#' Plot function for 'network' object
#' 
#' This function plots the given network. If the result of another plot is 
#' provided, this plot will be modified for easier comparison.
#' @param network A 'network' object.
#' @param compare_graph The plot of another network to use for comparison.
#' @param show_modules If TRUE, the modules will highlighted in the graph. 
#' Defaults to FALSE if there is exactly one module in the network and to TRUE
#' otherwise.
#' @param as_subgraph If TRUE, only nodes of positive degree will be shown. 
#' Defaults to FALSE if there are 100 or fewer nodes in the network and to TRUE
#' otherwise.
#' @param ... Additional arguments passed to plot_modules() or plot_network().
#' @return Creates a plot of the module and returns a graph object. 
#' See ?plot_modules and ?plot_network for details.
#' @return A 'network_plot' object for the network. This object can be passed 
#' back into a future call of plot.network() through the `compare_graph` 
#' argument, which will setup the plot for easier comparison between the old 
#' graph and the new graph of `network`.
#' @export
plot.network <- function(network, 
                         compare_graph = NULL, 
                         show_modules = FALSE, 
                         as_subgraph = FALSE, 
                         ...) {
  if(show_modules) {
    plot_modules(network, 
                 compare_graph = compare_graph, 
                 as_subgraph = as_subgraph, ...)
  } else {
    plot_network(network,
                 compare_graph = compare_graph, 
                 as_subgraph = as_subgraph, ...)
  }
}

#' Plot function for 'network_module' object.
#' 
#' @param module A 'network_module' object.
#' @param ... Additional arguments passed to plot_network().
#' @return Creates a plot of the module and returns a graph object. 
#' See ?plot_network for details.
#' @export
plot.network_module <- function(module, ...) {
  plot_network(module, ...)
}


#' Plot function for 'network_plot' class
#' 
#' @param network_plot A 'network_plot' object obtained from plot.network() or
#' plot_network().
#' @param ... Additional arguments passed to plot.igraph().
#' @export
plot.network_plot <- function(network_plot, 
                              vertex.color = network_plot$vertex.color, 
                              vertex.label.font = 2,
                              vertex.size = network_plot$vertex.size, 
                              vertex.frame.color = network_plot$vertex.frame.color,
                              vertex.label.color = network_plot$vertex.label.color, 
                              vertex.label.cex = 0.7,
                              edge.color = network_plot$edge.color, 
                              edge.width = network_plot$edge.width, 
                              layout = network_plot$coords,
                              ...) {
  
  plot(network_plot$graph, 
       vertex.color = vertex.color, 
       vertex.label.font = vertex.label.font,
       vertex.size = vertex.size, 
       vertex.frame.color = vertex.frame.color,
       vertex.label.color = vertex.label.color, 
       vertex.label.cex = vertex.label.cex,
       edge.color = edge.color, 
       edge.width = edge.width, 
       layout = layout,
       ...)
}


#' Plot matrix representation of network
#' 
#' This function plots the given network as a heatmap to visualize the 
#' corresponding association matrix. A weighted association matrix can be 
#' provided to show relative strength of associations.
#' @param network Either a network object or association matrix of the network.
#' @param main A string containing the title for the graph.
#' @param col Color palatte used for heatmap. See ?heatmap for details.
#' @return The matrix used to create the heatmap
#' @export
plot_network_matrix <- function(network, main = "Untitled",
                                col = colorRampPalette(RColorBrewer::brewer.pal(8, "Greys"))(50)) {
  if(class(network) == "network") {
    adj_matrix <- get_adjacency_matrix(network)
  } else if(is.matrix(network)){
    adj_matrix <- network
  } else {
    stop("network should be either a matrix or a network class object.")
  }
  
  heatmap(adj_matrix, main = main,
          symm = TRUE, Rowv = NA, Colv = NA, col = col)
  
  return(adj_matrix)
}






#' Plot the difference between two networks
#' 
#' This function plots the difference in connectivity between two networks. 
#' For two identical networks, the graph will be empty. For non-identical 
#' networks, black edges are shared by both networks but differ in magnitude or 
#' direction (if the networks are weighted), tan edges are in network_1 but not 
#' network_2, and red edges are in network_2 but not network_1. All edges are
#' scaled according to the strongest association in either network.
#' @param network_1 A 'network' or 'matrix' object.
#' @param network_2 A 'network' or 'matrix' object.
#' @param compare_graph The plot of another network to use for comparison.
#' @param as_subgraph If TRUE, only nodes of positive degree will be shown.
#' @param node_scale Used for scaling of nodes.
#' @param edge_scale Used for scaling of edges.
#' @param node_color The color used for the nodes.
#' @param edge_colors A vector of three colors used for edges; the first colors
#' edges common to both network, the second colors edges in network_1 but not
#' network_2, and the third colors edges that are in network_2 but not 
#' network_1. Default is c("black", "wheat", "red").
#' @param generate_coords A function to generate the layout of a graph; used
#' if coords is NULL. See ?igraph::layout_ for details. Other options include 
#' 'igraph::as_star', 'igraph::in_circle', and 'igraph::with_fr', among many others.
#' @param include_vertex_labels If TRUE, the verticies will be labeled.
#' @param display_plot If TRUE (default), the plot will be generated and displayed.
#' @param ... Additional arguments passed to plot.igraph().
#' @return Creates a plot of the network and returns a graph object. 
#' The graph object can be passed back into a future call of 'plot_network()',
#' 'plot_network_diff()' or 'plot_network_sim()'
#' through the 'compare_edge' argument, which will setup the plot for easier 
#' comparison between the old graph and the graph of 'network'.
#' @export
plot_network_diff <- function (network_1, network_2, compare_graph = NULL,
                               as_subgraph = FALSE,
                               node_scale = 4, edge_scale = 1, 
                               node_color = adjustcolor("orange", 0.5),
                               edge_colors = c("black", "wheat", "red"),
                               generate_layout = igraph::nicely,
                               include_vertex_labels = TRUE, 
                               ...) {
  ##################################
  # Check arguments for errors.
  
  if (!(class(network_1) %in% c("network", "network_module", "matrix")))
    stop(paste0("Argument 'network_1' must be a 'network', 'network_module', ",
                "or 'matrix' object."))
  
  if(!is.null(compare_graph)) {
    if(class(compare_graph) != "network_plot")
      stop("Argument 'compare_graph' must be a 'network_plot' object.")
    if((nrow(compare_graph$coords) != length(get_node_names(network_1))) &&
       !all(attr(igraph::V(compare_graph$graph), "names") %in% get_node_names(network_1)))
      stop(paste("Argument 'compare_graph' and 'network_1' must contain the same",
                 "number of nodes or contain an overlapping set of nodes."))
  }
  
  if (!(class(network_2) %in% c("network", "network_module", "matrix")))
    stop(paste0("Argument 'network_2' must be a 'network', 'network_module', ",
                "or 'matrix' object."))
  ##################################
  
  adj_matrix_1 <- get_adjacency_matrix(network_1)
  adj_matrix_2 <- get_adjacency_matrix(network_2)
  
  if(ncol(adj_matrix_1) != ncol(adj_matrix_2)) {
    stop("Arguments 'network_1' and 'network_2' must contain the same number of nodes.")
  }
  
  if(is.null(colnames(adj_matrix_1)) && is.null(colnames(adj_matrix_2))) {
    warning("Columns are unnamed in the adjacency matrix for 'network_1' and 'network_2")
    colnames(adj_matrix_1) <- 1:ncol(adj_matrix_1)
    colnames(adj_matrix_2) <- colnames(adj_matrix_1)
  } else if(is.null(colnames(adj_matrix_1))) {
    # If a matrix without column names is provded, we assume the columns
    # align with those from 'network_1'.
    warning(paste("Columns are unnamed in the adjacency matrix for 'network_1''",
                  "using the names from 'network_2'."))
    colnames(adj_matrix_1) <- colnames(adj_matrix_2)
  } else if(is.null(colnames(adj_matrix_2))) {
    warning(paste("Columns are unnamed in the adjacency matrix for 'network_2''",
                  "using the names from 'network_1'."))    
    colnames(adj_matrix_2) <- colnames(adj_matrix_1)
  } else if(!all(colnames(adj_matrix_1) == colnames(adj_matrix_2))) {
    stop("Arguments 'network_1' and 'network_2' must contain the same nodes (column names).")
  }
  
  weighted <- is_weighted(network_1) && is_weighted(network_2)
  # Set up association matricies if networks are weighted.
  if(weighted) {
    assoc_matrix_1 <- get_association_matrix(network_1)
    assoc_matrix_2 <- get_association_matrix(network_2)
    colnames(assoc_matrix_1) <- colnames(adj_matrix_1)
    colnames(assoc_matrix_2) <- colnames(adj_matrix_2)
    assoc_matrix_diff <- abs(assoc_matrix_1 - assoc_matrix_2)
    assoc_matrix_diff[assoc_matrix_diff < 10^-13] <- 0
    g <- igraph::graph_from_adjacency_matrix(assoc_matrix_diff, 
                                             mode = "undirected", 
                                             weighted = TRUE)
  } else {
    assoc_matrix_1 <- NULL
    assoc_matrix_2 <- NULL
    g <- igraph::graph_from_adjacency_matrix((adj_matrix_1 | adj_matrix_2) * 1, 
                                             mode = "undirected", 
                                             weighted = NULL)
  }
  
  # Initialize the comparison plot, if one is provided.
  if(is.null(compare_graph)) {
    g_compare <- NULL
  } else {
    # Initialize the graph to compare to. 
    g_compare <- compare_graph$graph
    # Subset both graphs to common nodes.
    common_nodes <- intersect(attr(igraph::V(g), "names"), attr(igraph::V(g_compare), "names"))
    index_subset_g <- which(attr(igraph::V(g), "names") %in% common_nodes)
    index_subset_g_compare <- which(attr(igraph::V(g_compare), "names") %in% common_nodes)
    
    # Update based on comparison graph.
    g <- igraph::induced_subgraph(g, index_subset_g)
    if(weighted) {
      assoc_matrix_1 <- assoc_matrix_1[index_subset_g, index_subset_g]
      assoc_matrix_2 <- assoc_matrix_2[index_subset_g, index_subset_g]
      assoc_matrix_diff <- assoc_matrix_diff[index_subset_g, index_subset_g]
    }
    adj_matrix_1 <- adj_matrix_1[index_subset_g, index_subset_g]
    adj_matrix_2 <- adj_matrix_2[index_subset_g, index_subset_g]
  }
  
  # Create subgraph, if requested.
  if(as_subgraph) {
    # Determine which nodes to keep in subgraph.
    index_nodes <- unname(which(igraph::degree(g) > 0))
    
    # If g contains no edges, then a subgraph cannot be created.
    if(length(index_nodes) == 0) {
      warning("Cannot create subgraph; the network must contain at least one edge.")
    } else {
      # Update based on subset of nodes.
      if(!is.null(g_compare)) {
        index_subset_g_compare <- index_subset_g_compare[index_nodes]
      }
      g <- igraph::induced_subgraph(g, index_nodes)
      if(weighted) {
        assoc_matrix_1 <- assoc_matrix_1[index_nodes, index_nodes]
        assoc_matrix_2 <- assoc_matrix_2[index_nodes, index_nodes]
        assoc_matrix_diff <- assoc_matrix_diff[index_nodes, index_nodes]
      }
      adj_matrix_1 <- adj_matrix_1[index_nodes, index_nodes]
      adj_matrix_2 <- adj_matrix_2[index_nodes, index_nodes]
    }
  }
  
  # Initialize coordinates for graph layout.
  if(is.null(compare_graph)) {
    coords = igraph::layout_(g, generate_layout())
  } else {
    coords = compare_graph$coords[index_subset_g_compare, ]
  }
  if(nrow(coords) != igraph::vcount(g)) {
    stop("coords do not match number of verticies in the graph.")
  }
  
  # Adjust node size, color, and frame.
  if(weighted) {
    # Scale associations relative to largest association in either network.
    temp <- max(c(abs(assoc_matrix_1[lower.tri(assoc_matrix_1)]),
                  abs(assoc_matrix_2[lower.tri(assoc_matrix_2)])))
    node_weights <- sqrt(apply(assoc_matrix_diff, 2, sum) /
                           ifelse(temp == 0, 1, temp))
  } else {
    # Default node weights are proportional to sqrt(degree change)
    node_weights <- sqrt(apply(abs(adj_matrix_1 - adj_matrix_2), 2, sum))
  }
  vertex.size <- node_weights * node_scale
  vertex.frame.color <- "white"
  if(include_vertex_labels) {
    vertex.label.color <- rgb(0.1, 0.1, 0.1, 0.8)
  } else {
    vertex.label.color <- rgb(0, 0, 0, 0) # Make labels transparent.
  }
  
  # Adjust edge width and color.
  if(weighted) {
    edge_weights <- assoc_matrix_diff[lower.tri(assoc_matrix_diff)]
    edge_weights <- edge_weights[edge_weights != 0]
    # Scale edges relative to largest association in either network.
    temp <- max(c(abs(assoc_matrix_1[lower.tri(assoc_matrix_1)]),
                  abs(assoc_matrix_2[lower.tri(assoc_matrix_2)])))
    edge_weights <- edge_weights / ifelse(temp == 0, 1, temp)
    if(length(edge_weights) != length(igraph::E(g))) {
      stop("Edge weights do not match number of edges.")
    }
  } else {
    edge_weights <- rep(1, length(igraph::E(g)))
  }
  
  edge.width <- edge_weights * edge_scale
  
  # Color is based on differences between graphs.
  g_1 <- igraph::graph_from_adjacency_matrix(adj_matrix_1, 
                                             mode = "undirected", 
                                             weighted = NULL)
  g_2 <- igraph::graph_from_adjacency_matrix(adj_matrix_2, 
                                             mode = "undirected", 
                                             weighted = NULL)
  edge.color <- vector("character", length(igraph::E(g)))
  # Edges in g1 are "wheat", in g2 are "red", and in both are "black".
  subset_1 <- which(attr(igraph::E(g), "vnames") 
                    %in% attr(igraph::E(g_1), "vnames"))
  edge.color[subset_1] <- edge_colors[2]
  subset_2 <- which(attr(igraph::E(g), "vnames") 
                    %in% attr(igraph::E(g_2), "vnames"))
  edge.color[subset_2] <- edge_colors[3]
  edge.color[intersect(subset_1, subset_2)] <- edge_colors[1]
  
  plot(g, vertex.color = node_color, vertex.label.font = 2,
       vertex.size = vertex.size, vertex.frame.color = vertex.frame.color,
       vertex.label.color = vertex.label.color, vertex.label.cex = 0.7,
       edge.color = edge.color, edge.width = edge.width, layout = coords, 
       ...)
  
  plot_summary <- list(graph = g,
                       coords = coords,
                       node_weights = node_weights,
                       edge_weights = edge_weights,
                       vertex.size = vertex.size,
                       vertex.frame.color = vertex.frame.color,
                       vertex.label.color = vertex.label.color,
                       vertex.color = node_color,
                       edge.color = edge.color,
                       edge.width = edge.width)
  
  class(plot_summary) <- "network_plot"
  
  # Return list silently:
  invisible(plot_summary)
}



#' Plot the similarity between two networks
#' 
#' This function plots the similarity of connections between two networks. 
#' Both networks must be weighted. The width of each edge corresponds to 
#' the strength of similarity and is calculated by sqrt(abs((s1 + s2)s1s2)), 
#' where s1 and s2 are the weights for a particular
#' connection in network_1 and network_2, respectively
#' @param network_1 A weighted 'network' or 'matrix' object.
#' @param network_2 A weighted 'network' or 'matrix' object.
#' @param compare_graph The plot of another network to use for comparison.
#' @param ... Additional arguments passed to 'plot_network()'.
#' @return Creates a plot of the network and returns a graph object. 
#' The graph object can be passed back into a future call of 'plot_network()',
#' 'plot_network_diff()' or 'plot_network_sim()'
#' through the 'compare_edge' argument, which will setup the plot for easier 
#' comparison between the old graph and the graph of 'network'.
#' @export
plot_network_sim <- function (network_1, network_2, compare_graph = NULL, ...) {
  ##################################
  # Check arguments for errors.
  
  if (!(class(network_1) %in% c("network", "network_module", "matrix")))
    stop(paste0("Argument 'network_1' must be a 'network', 'network_module', ",
                "or 'matrix' object."))
  
  if (!(class(network_2) %in% c("network", "network_module", "matrix")))
    stop(paste0("Argument 'network_2' must be a 'network', 'network_module', ",
                "or 'matrix' object."))
  ##################################
  
  if(!is_weighted(network_1) || !is_weighted(network_2)) {
    stop("Both networks must be weighted.")
  }
  
  A <- get_association_matrix(network_1)
  B <- get_association_matrix(network_2)
  
  if(ncol(A) != ncol(B)) {
    stop("Both networks must contain the same number of nodes.")
  }
  
  if(is.null(colnames(A)) && is.null(colnames(B))) {
    warning("Networks have unnamed nodes. Setting as 1:p.")
    colnames(A) <- 1:ncol(A)
    colnames(B) <- colnames(A)
  } else if(is.null(colnames(A))) {
    # If a matrix without column names is provded, we assume the columns
    # align with those from 'network_1'.
    warning(paste("Columns are unnamed 'network_1';",
                  "using the names from 'network_2'."))
    colnames(A) <- colnames(B)
  } else if(is.null(colnames(B))) {
    warning(paste("Columns are unnamed for 'network_2';",
                  "using the names from 'network_1'."))    
    colnames(B) <- colnames(A)
  } else if(!all(colnames(A) == colnames(B))) {
    stop("The networks must contain the same node names.")
  }
  
  genes <- colnames(A)
  p <- nrow(A)
  A <- A[lower.tri(A)]
  B <- B[lower.tri(B)]
  index <- which(abs(A) > 10^-2 & abs(B) > 10^-2)
  if(length(index) == 0) 
    return(matrix(0, p, p))
  S <- rep(0, length(A))
  A <- A[index]
  B <- B[index]
  S[index] <- sign(A + B) * sqrt(abs((A + B) * A * B))
  sim <- matrix(0, p, p)
  sim[lower.tri(sim)] <- S
  sim <- sim + t(sim)
  colnames(sim) <- genes
  
  plot_network(sim, compare_graph = compare_graph, ...)
}

#' Scatter plot of two gene expressions
#' 
#' Plots the expression of two genes for visual assessment of association.
#' @param counts_list A named list containing one or more n by p gene expression 
#' profiles, one for each group or subpopulation under consideration.
#' @param geneA The name of the first gene to plot. Must be either a character
#' string matching a column name in each matrix of counts_list or an integer
#' to index the columns.
#' @param geneB The name of the second gene to plot. Must be either a character
#' string matching a column name in each matrix of counts_list or an integer
#' to index the columns.
#' @param method Charater string either "lm" or "loess" used for plotting. 
#' For no line, set method = NULL. 
#' @param se_alpha Sets transparancy of confidence interval around association 
#' trend line. Set to 0 to remove the confidence interval.
#' @param do_facet_wrap If TRUE, the groups are plotted in seperate graphs.
#' @param scales Only used if do_facet_wrap is TRUE. See ggplot2::facet_wrap
#' for details.
plot_gene_pair <- function(counts_list, geneA, geneB, method = "loess", se_alpha = 0.1,
                           do_facet_wrap = FALSE, scales = "fixed") {
  if(!is.list(counts_list)) {
    counts_list <- list(x = counts_list)
    names(counts_list) <- "NA"
  } 
  
  groups <- names(counts_list)
  if(is.null(groups)) {
    groups <- as.character(1:length(counts_list))
  }
  exprA <- lapply(counts_list, function(x) {
    if(is.character(geneA)) {
      indexA <- which(colnames(x) == geneA)
    } else {
      indexA <- geneA
    }
    x[, indexA]
  })
  exprB <- lapply(counts_list, function(x) {
    if(is.character(geneB)) {
      indexB <- which(colnames(x) == geneB)
    } else {
      indexB <- geneB
    }
    x[, indexB]
  })
  n <- sapply(counts_list, nrow)
  
  df <- tibble::tibble(
    group = unlist(lapply(1:length(n), function(i) rep(groups[i], n[i]))),
    geneA = unlist(exprA),
    geneB = unlist(exprB))
  
  
  # Create the plot.
  g <- ggplot2::ggplot(df, ggplot2::aes(x = geneA, y = geneB, color = group))
  if(do_facet_wrap) {
    g <- g + ggplot2::facet_wrap(. ~ group, scales = scales)
  } 
  g <- g + ggplot2::geom_point(alpha = 0.5)
  if(!is.null(method)) {
    g <- g + ggplot2::geom_smooth(method =  method, alpha = se_alpha)
  }
  
  g <- g +
    ggplot2::theme_bw() +
    ggplot2::labs(x = paste("Expression of", geneA), 
                  y = paste("Expression of", geneB),
                  color = "Group")
  
  
  plot(g)
  
  return(g)
}
