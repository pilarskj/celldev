# Plot example trees

library(ape)
library(treeio)
library(ggtree)

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

tree_ids = c("tree_s_1", "tree_bd_15")
trueTrees = sapply(tree_ids, function(x) read.tree(paste0(dir, "/Trees/", x, ".newick")), simplify = F, USE.NAMES = T)

path = paste0(dir, "/Typewriter/baseline/inferenceOutput/")
mccTrees = sapply(tree_ids, function(x) read.beast(paste0(path, x, ".mcc.tree")), simplify = F, USE.NAMES = T)
posTrees = sapply(tree_ids, function(x) read.beast(paste0(path, x, ".thinned.trees")), simplify = F, USE.NAMES = T)

# origin/ experiment duration
Texp = 40

plot_trees <- function(x, title) {
  
  true = trueTrees[[x]]
  mcc = mccTrees[[x]]
  pos = posTrees[[x]]
  
  g1 = ggtree(true, root.position = Texp - max(node.depth.edgelength(true))) +
    geom_tiplab(size = 1) +
    expand_limits(x = c(0, 40))
  
  g2 = ggdensitree(pos, alpha = 0.01, tip.order = rev(get_taxa_name(g1))) +
    geom_tiplab(size = 1)
  
  g3 = ggtree(mcc, root.position = Texp - max(node.depth.edgelength(mcc@phylo))) +
    geom_range(range = 'height_0.95_HPD', center = 'height_median', color = 'lightblue', alpha = 0.5, size = 1.5) +
    geom_nodepoint(aes(color = posterior), size = 2) +
    geom_tiplab(size = 1) +
    expand_limits(x = c(0, 40)) +
    scale_color_continuous(name = "posterior\nsupport", limits = c(0,1), breaks = c(0, 0.5, 1), low = "darkred", high = "darkgreen") +
    theme_tree2() +
    theme(legend.position = c(0.1, 0.8), legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
          legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0.1, "cm"))
  
  g = (g1 + ggtitle(title) + theme(plot.title.position = "plot")) / g2 / g3 &
    theme(text = element_text(family = 'lmroman'))
  return(g)
}

# combined plot
g_syn = plot_trees(tree_ids[1], "A: synchronous tree")
g_bd = plot_trees(tree_ids[2], "B: birth-death tree with sampling")

g_syn | g_bd
ggsave("pdf/SuppFigure2_homTrees.pdf", height = 11, width = 8)
