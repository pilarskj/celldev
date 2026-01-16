# -----
# Script for tree simulations (for homogeneous cell populations) - parameters adapted for generating larger trees
# -----

library(ape)
library(TreeSim)
library(stringr)
library(dplyr)
library(treeio) # for loading .trees file

# set path to repo as working directory
setwd("~/Projects/celldev")

# load functions
source("homog/simulation_trees/tree_functions.R")

# specify output paths
params_f = 'homog/simulation_trees/large_tree_params.csv'
tree_dir = '~/Projects/celldev_data/homog/LargeTrees'
if (!dir.exists(tree_dir)) {dir.create(tree_dir)}

# specify parameters
Texp = 40 # experiment duration (tree origin)
Nsim = 20 # number of simulations

treeParams = list()


# SIMULATE TREES

# complete synchronous tree 
tree_s = generateSynTree(divisions = 10, Texp = Texp, rho = 1)
plot(tree_s)
axisPhylo()

# store tree parameters
treeParams[[1]] = data.frame("tree" = paste0("tree_s_", c(1:Nsim)), "nDiv" = 10, "deathProb" = 0, "rho" = 1, 
                             "nCells" = rep(length(tree_s$tip.label), Nsim), 
                             "treeLength" = rep(tree_s$length, Nsim), "treeHeight" = rep(tree_s$height, Nsim))


# synchronous trees with sampling
trees_ss = lapply(c(1:Nsim), function(i) {generateSynTree(divisions = 13, Texp = Texp, rho = 0.1, seed = i)})
sapply(trees_ss, Ntip) %>% unique # 819

# store tree parameters
treeParams[[2]] = data.frame("tree" = paste0("tree_ss_", c(1:Nsim)), "nDiv" = 13, "deathProb" = 0, "rho" = 0.1,
                             "nCells" = sapply(trees_ss, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_ss, function(t) {t$length}), 
                             "treeHeight" = sapply(trees_ss, function(t) {t$height}))


# synchronous trees with cell death
trees_sd = list()
i = 1
while (length(trees_sd) < Nsim) {
  # simulate tree
  sd <- generateSynDeathTree(deathprob = 0.15, divisions = 13, Texp = Texp, rho = 1, seed = i)
  # print(Ntip(sd))
  if (length(sd$tip.label) %in% c(200:2000)) {
    # save if number of tips within range
    trees_sd = append(trees_sd, list(sd))
  } 
  i = i + 1
}
print(i) # number of simulations needed: 22
sapply(trees_sd, Ntip) %>% summary 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 367.0   588.5   971.5  1039.2  1499.5  1668.0 

# store tree parameters
treeParams[[3]] = data.frame("tree" = paste0("tree_sd_", c(1:Nsim)), "nDiv" = 13, "deathProb" = 0.15, "rho" = 1,
                             "nCells" = sapply(trees_sd, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_sd, function(t) {t$length}), 
                             "treeHeight" = sapply(trees_sd, function(t) {t$height}))


# synchronous trees with cell death and sampling
trees_sds = list()
i = 1
while (length(trees_sds) < Nsim) {
  # simulate tree
  sds <- generateSynDeathTree(deathprob = 0.12, divisions = 16, Texp = Texp, rho = 0.1, seed = i)
  print(Ntip(sds))
  if (length(sds$tip.label) %in% c(200:2000)) {
    # save if number of tips within range
    trees_sds = append(trees_sds, list(sds))
  } 
  i = i + 1
}
print(i) # number of simulations needed: 21
sapply(trees_sds, Ntip) %>% summary 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 486.0   802.8   994.0   952.0  1098.2  1317.0 

# store tree parameters
treeParams[[4]] = data.frame("tree" = paste0("tree_sds_", c(1:Nsim)), "nDiv" = 16, "deathProb" = 0.12, "rho" = 0.1,
                             "nCells" = sapply(trees_sds, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_sds, function(t) {t$length}),
                             "treeHeight" = sapply(trees_sds, function(t) {t$height}))


# birth-death trees
# because of runtime, simulated with ReMASTER (using script remaster.xml)
# load .trees file
trees_bd = read.beast(file.path(tree_dir, "trees_bd.trees"))
sapply(trees_bd, Ntip) %>% summary
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 332.0   559.2   797.0   882.1  1229.5  1463.0 

# re-assign tip labels, calculate tree heights and lengths
trees_bd = lapply(trees_bd, function(tree) {
  tree <- tree@phylo
  tree$tip.label <- c(0:(length(tree$tip.label)-1))
  tree$length <- sum(tree$edge.length)
  tree$height <- max(node.depth.edgelength(tree))
  return(tree)
})

# store tree parameters
treeParams[[5]] = data.frame("tree" = paste0("tree_bd_", c(1:Nsim)), "birthRate" = 0.27, "deathRate" = 0.035, "rho" = 0.1,
                             "nCells" = sapply(trees_bd, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_bd, function(t) {t$length}), 
                             "treeHeight" = sapply(trees_bd, function(t) {t$height}))


# save trees
lapply(c(1:20), function(i) {
  write.tree(tree_s, file = file.path(tree_dir, paste0("tree_s_", i, ".newick")))
  write.tree(trees_ss[[i]], file = file.path(tree_dir, paste0("tree_ss_", i, ".newick")))
  write.tree(trees_sd[[i]], file = file.path(tree_dir, paste0("tree_sd_", i, ".newick")))
  write.tree(trees_sds[[i]], file = file.path(tree_dir, paste0("tree_sds_", i, ".newick")))
  write.tree(trees_bd[[i]], file = file.path(tree_dir, paste0("tree_bd_", i, ".newick")))
})


# approximation of birth and death rates for synchronous trees
treeParams[1:4] = lapply(treeParams[1:4], function(p) {
  # population-level 
  p %>% mutate(birthRate = (1-deathProb) * (nDiv / Texp) * log(2),
               deathRate = deathProb * (nDiv / Texp) * log(2))
  # per-lineage
  # p %>% mutate(birthRate = (1-deathProb) * ((nDiv+1) / Texp),
  #              deathRate = deathProb * ((nDiv+1) / Texp))
})

# calculate death probability for birth-death trees
treeParams[[5]] = treeParams[[5]] %>% mutate(deathProb = deathRate / (birthRate + deathRate))
treeParams = bind_rows(treeParams)

# save tree parameters
write.csv(treeParams, file = params_f, quote = F, row.names = F)

# # view stats
# treeParams = read.csv("~/Projects/celldev/homog/simulation_trees/large_tree_params.csv")
# treeParams %>% 
#   mutate(treeModel = str_extract(tree, "tree_[a-z]+")) %>%
#   group_by(treeModel) %>% 
#   summarise(min(nCells), max(nCells), median(nCells), mean(nCells), sd(nCells))


# --------------
# Sample cells to get smaller trees at lower sampling proportion (take another 10% of cells, yielding rho = 0.01 = 1%)
prop = 0.1
sampledTreeParams = list()

sample_tree <- function(tree) {
  sampled_tips = sample(tree$tip.label, round(prop*Ntip(tree), 0))
  tree = keep.tip(tree, tip = sampled_tips)
  tree$tip.label = c(0:(length(tree$tip.label)-1))
  tree$length = sum(tree$edge.length)
  tree$height = max(node.depth.edgelength(tree))
  return(tree)
}

set.seed(1)
sampled_trees_ss = lapply(trees_ss, sample_tree) %>% setNames(paste0("tree_ss_", c(1:Nsim)))
sampled_trees_sds = lapply(trees_sds, sample_tree) %>% setNames(paste0("tree_sds_", c(1:Nsim)))
sampled_trees_bd = lapply(trees_bd, sample_tree) %>% setNames(paste0("tree_bd_", c(1:Nsim)))
sampledTrees = c(sampled_trees_ss, sampled_trees_sds, sampled_trees_bd)

sampledTreeParams = lapply(sampledTrees, function(t) {
  df = data.frame("nCells" = length(t$tip.label),
                  "treeLength" = t$length, 
                  "treeHeight" = t$height) }) %>%
  bind_rows(.id = "tree") %>%
  left_join(treeParams %>% select(tree, nDiv, deathProb, birthRate, deathRate), by = "tree") %>%
  mutate(rho = 0.01) %>%
  select(all_of(colnames(treeParams))) # reorder columns

summary(sampledTreeParams$nCells)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 33.00   80.00   82.00   88.38  103.25  145.00 

# compare to previous trees (on average somewhat smaller)
test = read.csv("~/Projects/celldev/homog/simulation_trees/tree_params.csv")
test = test %>% 
  mutate(treeModel = str_extract(tree, "tree_[a-z]+")) %>%
  filter(treeModel %in% c("tree_ss", "tree_sds", "tree_bd"))
summary(test$nCells)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 25.00   82.00  102.00   99.28  123.25  180.00 

# save trees and parameters
sampled_params_f = 'homog/simulation_trees/sampled_tree_params.csv'
sampled_tree_dir = '~/Projects/celldev_data/homog/SampledTrees'
if (!dir.exists(tree_dir)) {dir.create(tree_dir)}

lapply(c(1:length(sampledTrees)), function(i) {
  write.tree(sampledTrees[[i]], file = file.path(sampled_tree_dir, paste0(names(sampledTrees)[i], ".newick")))
})

write.csv(sampledTreeParams, file = sampled_params_f, quote = F, row.names = F)
