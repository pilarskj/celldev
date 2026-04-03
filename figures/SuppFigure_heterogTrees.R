# Plot example multi-type trees

library(dplyr)
library(tidyr)
library(ape)
library(treeio)
library(ggtree)
library(ggplot2)
library(patchwork)
library(showtext)

setwd("~/Projects/celldev/figures")
dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-CellDev/heterog"

# plotting settings
palette = c("0" = "#56B4E9", "1" = "#009E73", "2" = "#E69F00", "3" = "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()

# trees
#tree_ids = c("distinct" = "tree_31", "hierarchical" = "tree_11")
tree_ids = c("tree_8", "tree_13", "tree_16", "tree_17")
names(tree_ids) = rep("hierarchical", 4)

trueTrees = sapply(c(1:length(tree_ids)), function(i) {
  f = paste0(dir, "/Trees/", names(tree_ids)[i], "/", tree_ids[i], ".newick")
  newick = readLines(f)
  tree = read.beast.newick(textConnection(paste0(newick, ";"))) 
  tree@data$type = as.factor(tree@data$type)
  return(tree)
}, simplify = F, USE.NAMES = T)

mccTrees = sapply(c(1:length(tree_ids)), function(i) {
  tree = read.beast(paste0(dir, "/Typewriter/", names(tree_ids)[i], "/inferenceOutput/", tree_ids[i], ".mcc.tree"))
  tree@data$type = as.factor(tree@data$type)
  return(tree)
  }, simplify = F, USE.NAMES = T)

# posTrees = sapply(c(1:length(tree_ids)), function(i) {
#   trees = read.beast(paste0(dir, "/Typewriter/", names(tree_ids)[i], "/inferenceOutput/", tree_ids[i], ".thinned.trees"))
#   trees = lapply(trees, function(tree) {
#     tree@data$type = as.factor(tree@data$type)
#     return(tree)
#   })
#   return(trees)
#   }, simplify = F, USE.NAMES = T)


# origin/ experiment duration
Texp = 40


# plot 
plot_trees <- function(x, title = NULL) {
  
  true = trueTrees[[x]]
  mcc = mccTrees[[x]]
  #pos = posTrees[[x]]
  
  g1 = ggtree(true, aes(color = type), root.position = Texp - max(node.depth.edgelength(true@phylo))) +
    geom_tiplab(size = 1) +
    expand_limits(x = c(0, 40)) +
    scale_color_manual(name = "Type", values = palette)
  
  # g2 = ggdensitree(pos, aes(color = type), alpha = 0.01, tip.order = rev(get_taxa_name(g1))) +
  #   geom_tiplab(size = 1) +
  #   scale_color_manual(values = palette) +
  #   theme(legend.position = "none")
  
  # rotate tree
  mcc@phylo = rotateConstr(mcc@phylo, rev(get_taxa_name(g1)))
  
  # pie charts on each node
  pie_data = mcc %>%
    as_tibble() %>%
    filter(is.na(label)) %>% # keep only internal nodes
    select(node, type.set, type.set.prob) %>%
    unnest(c(type.set, type.set.prob)) %>%
    mutate(type.set = as.factor(type.set)) %>%
    pivot_wider(names_from = type.set, values_from = type.set.prob, values_fill = 0)
  pies = nodepie(pie_data, cols = c("0", "1", "2", "3"))
  pies = lapply(pies, function(g) g + scale_fill_manual(values = palette))

  g3 = ggtree(mcc, root.position = Texp - max(node.depth.edgelength(mcc@phylo))) +
    #geom_range(range = 'height_0.95_HPD', center = 'height_median', color = 'lightgrey', alpha = 0.5, size = 1.5) +
    geom_inset(pies, width = 0.05, height = 0.05) +
    geom_tiplab(aes(color = type), size = 1) +
    expand_limits(x = c(0, 40)) +
    scale_color_manual(values = palette) +
    #theme_tree2() +
    theme(legend.position = "none")
  
  g = g1 + g3 &
  #g = (g1 + ggtitle(title) + theme(plot.title.position = "plot")) / g2 / g3 &
     theme(text = element_text(family = 'lmroman'))
  return(g)
}

# combined plot
# g_d = plot_trees(1, "A: terminal transitions")
# g_h = plot_trees(2, "B: chain-like transitions")
# 
# (g_d | g_h) + plot_layout(guides = "collect") 
# ggsave("pdf/SuppFigure14_heterogTrees.pdf", height = 11, width = 9)

g = lapply(c(1:4), plot_trees)
wrap_plots(g, nrow = 4, guides = "collect")
ggsave("pdf/SuppFigure15_biasedTrees.pdf", height = 9, width = 8)




