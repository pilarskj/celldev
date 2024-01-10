# -----
# Functions for simulating synchronous trees with cell death 
# -----

dtree <- function(deathprob, divisions, T0) {
  "Generate synchronous tree with cell death"
  # generate complete tree
  n_tips <- 2**divisions
  length_branch <- T0/(divisions + 1)
  tree_complete <- stree(n_tips, type = "balanced", tip.label = c(1:n_tips))
  tree_complete$edge.length <- rep(length_branch, nrow(tree_complete$edge))
  
  # extract edge matrix, root and inner nodes
  edge <- tree_complete$edge
  edge.length <- tree_complete$edge.length
  
  # go along branches, let cells die with probability deathprob
  i <- 1
  while (i <= nrow(edge)) {
    branch <- edge[i, ]
    die <- sample(c(TRUE, FALSE), 1, prob = c(deathprob, 1 - deathprob))
    if (die) {
      # shorten branch to time of cell death 
      edge.length[i] <- runif(1, min = 0, max = length_branch)
      dying_node <- branch[2]
      if (!dying_node %in% tree_complete$tip.label) {
        # remove all offspring of dying cell 
        ix_prev <- i + which(edge[-c(1:i), 1] < dying_node)[1]
        if (!is.na(ix_prev)) {
          edge <- edge[-c((i + 1) : (ix_prev - 1)), ] 
          edge.length <- edge.length[-c((i + 1) : (ix_prev - 1))] 
        } else {
          edge <- edge[-c((i + 1) : nrow(edge)), ] 
          edge.length <- edge.length[-c((i + 1) : length(edge.length))] 
        }
      }
    }
    i <- i + 1
  }
  # restore phylo object 
  tree <- restoreTree(edge, edge.length)
  return(tree)
  # for plotting surviving cells in color over whole tree
  #ix_in <- which(apply(tree_balanced$edge, 1, function(x) all(x %in% edge)))
  #edge_color <- rep("black", nrow(tree_balanced$edge))
  #edge_color[ix_in] <- "blue"
  #plot(tree_balanced, edge.color = edge_color)
  #axisPhylo()
}


restoreTree <- function(edge, edge.length) {
  "Helper function for restoring phylo object after cell death"
  # find and rename tips
  ix_tips <- which(!(edge[, 2] %in% edge[, 1]))
  tip.label <- c(1:length(ix_tips))
  edge[ix_tips, 2] <- tip.label
  
  # find and rename inner nodes
  innodes <- unique(edge[, 1])
  Nnode <- length(innodes)
  new_innodes <- c((length(ix_tips) + 1):(length(ix_tips) + length(innodes)))
  innodes <- setNames(innodes, new_innodes)
  edge[, 1] <- sapply(edge[, 1], function(x) {x <- as.integer(names(which(innodes == x))) }, simplify = T)
  if (!is.null(edge[-ix_tips, ]) & nrow(edge[-ix_tips, ]) != 0) {
    edge[-ix_tips, 2] <- sapply(edge[-ix_tips, 2], function(x) {x <- as.integer(names(which(innodes == x))) }, simplify = T)
  }
  
  # define updated phylo object
  phy <- list(edge = edge, edge.length = edge.length, tip.label = tip.label, Nnode = Nnode)
  class(phy) <- "phylo"
  attr(phy, "order") <- "cladewise"
  return(phy)
}

