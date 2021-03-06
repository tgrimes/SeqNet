% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/networks.R
\name{create_modules_for_network}
\alias{create_modules_for_network}
\title{Randomly sample subsets of genes for each module}
\usage{
create_modules_for_network(n_modules, p, avg_module_size = 50,
  sd_module_size = 50, min_module_size = 10, max_module_size = 200,
  sample_link_nodes_fn = sample_link_nodes,
  sample_module_nodes_fn = sample_module_nodes, ...)
}
\arguments{
\item{n_modules}{The number of modules to include in the network.}

\item{p}{The number of nodes in the network.}

\item{avg_module_size}{The average number of nodes in a module.}

\item{sd_module_size}{The standard deviation of module size.}

\item{min_module_size}{The minimum number of nodes in a module.}

\item{max_module_size}{A positive value. Any generated module sizes above this 
value will be reduced to 'max_module_size'. Set to 'Inf' to avoid this 
truncation.}

\item{sample_link_nodes_fn}{A function used for sampling link nodes for a new 
module.}

\item{sample_module_nodes_fn}{A function used for sampling nodes for a new 
module.}

\item{...}{Additional arguments passed to random_module().}
}
\value{
A list containing the indicies for genes contained in each module.
}
\description{
Creates a collection of modules containing randomly samples genes.
}
\examples{
# Create a two modules (having random structures and sizes) from a pool 
# of 100 nodes.
create_modules_for_network(n_modules = 2, p = 100)
# Set n_modules = NULL to continue making modules until all nodes have
# been selected at least once.
create_modules_for_network(n_modules = NULL, p = 100)
}
