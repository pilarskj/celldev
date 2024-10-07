# -----
# Script for tree simulations (for homogeneous cell populations)
# -----

library(ape)
library(TreeSim)
library(stringr)
library(dplyr)

# set path to repo as working directory
# setwd("~/Projects/celldev")

# load functions
source("homog/simulation_trees/tree_functions.R")

# specify output paths
params_f = 'homog/simulation_trees/tree_params.csv'
tree_dir = '~/Projects/celldev_data/homog/Trees'
if (!dir.exists(tree_dir)) {dir.create(tree_dir)}

# specify parameters
Texp = 40 # experiment duration (tree origin)
Nsim = 20 # number of simulations

treeParams = list()


# SIMULATE TREES

# complete synchronous tree 
tree_s = generateSynTree(divisions = 7, Texp = Texp, rho = 1)
plot(tree_s)
axisPhylo()

# store tree parameters
treeParams[[1]] = data.frame("tree" = "tree_s", "nDiv" = 7, "deathProb" = 0, "rho" = 1, 
                             "nCells" = length(tree_s$tip.label), 
                             "treeLength" = tree_s$length, "treeHeight" = tree_s$height)


# synchronous trees with sampling
trees_ss = lapply(c(1:Nsim), function(i) {generateSynTree(divisions = 10, Texp = Texp, rho = 0.1, seed = i)})
plot(trees_ss[[1]]) # plot example
axisPhylo()

# store tree parameters
treeParams[[2]] = data.frame("tree" = paste0("tree_ss_", c(1:Nsim)), "nDiv" = 10, "deathProb" = 0, "rho" = 0.1,
                             "nCells" = sapply(trees_ss, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_ss, function(t) {t$length}), 
                             "treeHeight" = sapply(trees_ss, function(t) {t$height}))


# synchronous trees with cell death
trees_sd = list()
i = 1
while (length(trees_sd) < Nsim) {
  # simulate tree
  sd <- generateSynDeathTree(deathprob = 0.2, divisions = 10, Texp = Texp, rho = 1, seed = i)
  if (length(sd$tip.label) %in% c(20:200)) {
    # save if number of tips within [20,200]
    trees_sd = append(trees_sd, list(sd))
  } 
  i = i + 1
}
print(i) # number of simulations needed
plot(trees_sd[[1]]) # plot example
axisPhylo()

# store tree parameters
treeParams[[3]] = data.frame("tree" = paste0("tree_sd_", c(1:Nsim)), "nDiv" = 10, "deathProb" = 0.2, "rho" = 1,
                             "nCells" = sapply(trees_sd, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_sd, function(t) {t$length}), 
                             "treeHeight" = sapply(trees_sd, function(t) {t$height}))


# synchronous trees with cell death and sampling
trees_sds = list()
i = 1
while (length(trees_sds) < Nsim) {
  # simulate tree
  sds <- generateSynDeathTree(deathprob = 0.1, divisions = 12, Texp = Texp, rho = 0.1, seed = i)
  if (length(sds$tip.label) %in% c(20:200)) {
    # save if number of tips within [20,200]
    trees_sds = append(trees_sds, list(sds))
  } 
  i = i + 1
}
print(i) # number of simulations needed
plot(trees_sds[[1]]) # plot example
axisPhylo()

# store tree parameters
treeParams[[4]] = data.frame("tree" = paste0("tree_sds_", c(1:Nsim)), "nDiv" = 12, "deathProb" = 0.1, "rho" = 0.1,
                             "nCells" = sapply(trees_sds, function(t) {length(t$tip.label)}),
                             "treeLength" = sapply(trees_sds, function(t) {t$length}),
                             "treeHeight" = sapply(trees_sds, function(t) {t$height}))


# birth-death trees
trees_bd = list()
i = 1
while (length(trees_bd) < 20) {
  # simulate birth-death tree
  bd = generateBirthDeathTree(birth = 0.2, death = 0.025, Texp = Texp, rho = 0.1, seed = i)
  if (typeof(bd) == "list") {
    if (length(bd$tip.label) %in% c(20:200)) {
      # save if number of tips within [20,200]
      trees_bd = append(trees_bd, list(bd))
    }
  } 
  i = i + 1
}
print(i) # number of simulations needed
plot(trees_bd[[1]]) # plot example
axisPhylo()

# store tree parameters
treeParams[[5]] = data.frame("tree" = paste0("tree_bd_", c(1:Nsim)), "birthRate" = 0.2, "deathRate" = 0.025, "rho" = 0.1,
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


# calculate missing tree parameters
treeParams = bind_rows(treeParams)

# calculate approximated population-level birth and death rates for synchronous trees
treeParams[1:61, ] = treeParams[1:61, ] %>% mutate(birthRate = (1-deathProb) * (nDiv / Texp) * log(2),
                                                   deathRate = deathProb * (nDiv / Texp) * log(2))
# per-lineage approximation
# treeParams[1:61, ] = treeParams[1:61, ] %>% mutate(birthRate = (1-deathProb) * ((nDiv+1) / Texp),
#                                                    deathRate = deathProb * ((nDiv+1) / Texp))

# calculate death probability for birth-death trees
treeParams[62:81, ] = treeParams[62:81, ] %>% mutate(deathProb = deathRate / (birthRate + deathRate))

# save tree parameters
write.csv(treeParams, file = params_f, quote = F, row.names = F)