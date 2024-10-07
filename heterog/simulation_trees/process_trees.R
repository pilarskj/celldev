# -----
# Script for processing simulated BDMM trees 
# (selecting valid trees, extracting cell types, summarizing tree parameters)
# -----

library(dplyr)
library(stringr)
library(ape)

# set path to repo as working directory
# setwd("~/Projects/celldev")

tree_dir = '~/Projects/celldev_data/Trees'
if (!dir.exists(tree_dir)) {dir.create(tree_dir)}

# specify prototype and number of simulations performed
prototype = "hierarchical" # or "distinct" 
nSim = 20

# specify path to output file
params_f = paste0('heterog/simulation_trees/tree_params_', prototype, '.csv')

# initialize
v = 0
treeParams = data.frame(matrix(data = NA, nrow = nSim, ncol = 4))
colnames(treeParams) = c("seed", "ntips", "treeHeight", "treeLength")

for (i in c(1:nSim)) {
  
  # get tree in newick format
  newick = readLines(file.path(tree_dir, paste0(prototype, "/tree_", i, ".newick")))
  
  # read tree as phylo object
  tree = read.tree(text = paste0(str_remove_all(newick, "\\[&type=...\\]"), ";"))
  
  # get tips and their types
  tips = tree$tip.label
  tip_types = str_extract(newick, paste0("[^0-9]", tips, "\\[&type=...\\]:"))
  tip_types = str_extract(str_extract(tip_types, "=\"[0123]"), "[0123]")
  taxa = data.frame(tips, tip_types)
  
  # select valid trees
  if (prototype == "hierarchical") {
    # check if min. 1 tip is of type 3
    if ("3" %in% tip_types) {
      valid = TRUE
      v = v + 1
    } else {
      valid = FALSE
    }
  } else if (prototype == "distinct") {
    # check if min. 1 tip of each type is present
    if (all(c("1", "2", "3") %in% tip_types)) {
      valid = TRUE
      v = v + 1
    } else {
      valid = FALSE
    }
  }
  
  # exclude if more than 200 tips
  if (length(tips) > 200) {
    valid = FALSE
    v = v - 1
  }
  
  if (valid) {
    # write csv with types of taxa
    write.table(taxa, file = file.path(tree_dir, paste0(prototype, "/tree_", i, ".types.csv")), 
                sep = ",", row.names = F, col.names = F, quote = F)
    
    # store tree parameters
    treeParams[i, "seed"] = i
    treeParams[i, "ntips"] = length(tips)
    
    log = readLines(file.path(tree_dir, paste0(prototype, "/tree_", i, ".log")))[2] # read tree log
    treeParams[i, c("treeHeight", "treeLength")] = str_split(log, pattern = "\t", simplify = T)[, c(2,3)]
    
  } else {
    print(paste0("Tree ", i, " not valid"))
    # remove tree
    unlink(file.path(tree_dir, paste0(prototype, "/tree_", i, ".newick")))
    unlink(file.path(tree_dir, paste0(prototype, "/tree_", i, ".log")))
  }
}

print(paste0(v, " trees valid!"))

# save tree parameters
treeParams = treeParams %>% drop_na()
write.table(treeParams, file = params_f, 
            sep = ",", row.names = F, col.names = T, quote = F)
