# -----
# Script for processing simulated BDMM trees 
# (selecting valid trees, extracting cell types, summarizing tree parameters)
# -----

library(dplyr)
library(stringr)
library(ape)

# setwd("~/Desktop/Part2")

# specify prototype and number of simulations performed
prototype = "hierarchical"
nSim = 25

# initialize
v = 0
treeParams = data.frame(matrix(data = NA, nrow = nSim, ncol = 4))
colnames(treeParams) = c("seed", "ntips", "treeHeight", "treeLength")

for (i in c(1:nSim)) {
  
  # get tree in newick format
  newick = readLines(paste0("Trees/", prototype, "/tree_", i, ".newick"))
  
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
  
  if (valid) {
    # write csv with types of taxa
    write.table(taxa, file = paste0("Trees/", prototype, "/tree_", i, ".types.csv"), 
                sep = ",", row.names = F, col.names = F, quote = F)
    
    # store tree parameters
    treeParams[i, "seed"] = i
    treeParams[i, "ntips"] = length(tips)
    
    log = readLines(paste0("Trees/", prototype, "/tree_", i, ".log"))[2] # read tree log
    treeParams[i, c("treeHeight", "treeLength")] = str_split(log, pattern = "\t", simplify = T)[, c(2,3)]
    
  } else {
    print(paste0("Tree ", i, " not valid"))
    # remove tree
    unlink(paste0("Trees/", prototype, "/tree_", i, ".newick"))
    unlink(paste0("Trees/", prototype, "/tree_", i, ".log"))
  }
}

print(paste0(v, " trees valid!"))

# save tree parameters
treeParams = treeParams %>% drop_na()
write.table(treeParams, file = paste0("Trees/", prototype, "/treeParams.csv"), 
            sep = ",", row.names = F, col.names = T, quote = F)
