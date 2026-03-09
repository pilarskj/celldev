# ----------
# Script for the analysis of typed trees inferred using BDMM-Prime in BEAST2
# ----------

# load packages
library(ape)
library(treeio)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(HDInterval)

rm(list = ls()) # clear environment

# adapt paths
codeDir = '~/Projects/celldev'
dataDir = '/Projects/celldev_data/heterog'
setwd(codeDir)


### Function for collecting type transition statistics from trees
get_trans_stats <- function(tree, origin = 40, origin_type = "0", add_stem = F, count_tips = F) {
  # convert branching times
  height = max(node.depth.edgelength(tree@phylo))
  times = abs(node.depth.edgelength(tree@phylo) - height)
  df = tree %>% as_tibble() %>% mutate(time = times, type = as.character(type))
  df$isTip = isTip(tree, df$node)
  
  # optionally, add stem branch
  if (add_stem) {
    root = Ntip(tree) + 1
    origin_node = Ntip(tree) + Nnode(tree) + 1
    df[df$node == root, "parent"] = origin_node
    df[df$node == root, "branch.length"] = origin - height
    df = rbind(df, data.frame(parent = origin_node, node = origin_node, branch.length = 0, label = NA, type = origin_type, time = origin, isTip = F))
  }
  
  # find type transitions along branches
  ix = df %>% group_by(parent) %>% summarise(children = n()) %>% filter(children == 1) %>% pull(parent)
  trans = df %>% filter(parent %in% ix) %>% select(parent, node, type) %>%
    left_join(df %>% filter(node %in% ix) %>% select(parent = node, parent_type = type, time), by = "parent") %>%
    select(parent, node, parent_type, type, time)
  
  # compute time spent in each type
  stats = df %>% group_by(type) %>% summarise(time_spent = sum(branch.length)) #%>%
    # mutate(prop_spent = time_spent / sum(time_spent))
  # get number of tips per type
  if (count_tips) {
    stats = stats %>% full_join(df %>% filter(isTip) %>% group_by(type) %>% summarise(n_tips = n()), by = "type") #%>%
      # mutate(prop_tips = n_tips / sum(n_tips, na.rm = T))
  }
  # compute number of transitions to each type and identify timing of first transition
  stats = stats %>% full_join(trans %>% group_by(type) %>% summarise(n_trans = n(), first_trans = max(time)), by = "type")
  return(stats)
}


### TRUE TREES
# multi-type model prototype
prototype = "hierarchical"

# tree parameters
treeParams = read.csv(paste0(codeDir, "/heterog/simulation_trees/tree_params_", prototype, ".csv"))  
seeds = treeParams$seed

# read true trees
trueTrees = lapply(seeds, function(i) {
  f = paste0(dataDir, "/Trees/", prototype, "/tree_", i, ".newick")
  newick = readLines(f)
  tree = read.beast.newick(textConnection(paste0(newick, ";"))) 
  return(tree)
})
names(trueTrees) = seeds

# compute type transition statistics and save them
stats = lapply(trueTrees, get_trans_stats, add_stem = T, count_tips = T) %>% bind_rows(.id = "tree")
stats[is.na(stats$n_tips), "n_tips"] = 0
write.csv(stats, file = paste0(dataDir, "/Trees/", prototype, "/transition_stats.csv"), row.names = F, quote = F)


### POSTERIOR TREE STATISTICS
# lineage recorder
method = "TiDe"

# directory for current setup
dir = paste0(dataDir, "/", method, "/", prototype)

for (seed in seeds) {
  # read posterior trees
  postTrees = read.beast(paste0(dir, "/inferenceOutput/tree_", seed, ".typed.trees"))
  # remove 10% burnin
  burnin = ceiling(length(postTrees) * 0.1)
  postTrees = postTrees[c((burnin+1):length(postTrees))]
  # calculate tree statistics
  postStats = lapply(postTrees, get_trans_stats) %>% bind_rows(.id = "Sample")
  # store as .log file
  write.table(postStats, file = paste0(dir, "/inferenceOutput/tree_", seed, ".typestats.log"), 
              quote = F, sep = "\t", row.names = F)
  print(paste0("tree ", seed, " done"))
}


### SUMMARY TREE STATISTICS
hpd_lower <- function(x) hdi(na.omit(x), credMass = 0.95)[1]
hpd_upper <- function(x) hdi(na.omit(x), credMass = 0.95)[2]

treeStats = lapply(seeds, function(seed) {
  postStats = read.table(file = paste0(dir, "/inferenceOutput/tree_", seed, ".typestats.log"), 
                         header = T, sep = "\t")
  
  sumStats = postStats %>%
    group_by(type) %>%
    summarise(across(
      .cols = c(time_spent, n_trans, first_trans),
      .fns = list(median = ~median(.x, na.rm = T), lower = hpd_lower, upper = hpd_upper),
      .names = "{.col}_{.fn}"))
  
  return(sumStats)
})
names(treeStats) = seeds
treeStats = treeStats %>% bind_rows(.id = "tree")

saveRDS(treeStats, file = paste0(dir, "/analysisOutput/treeStats.Rdat"))

