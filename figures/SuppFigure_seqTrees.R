# Plot example trees

library(ape)
library(treeio)
library(ggtree)

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

tree_ids = c("tree_s_1", "tree_ss_1", "tree_sd_2", "tree_sds_8", "tree_bd_15")
trueTrees = sapply(tree_ids, function(x) read.tree(paste0(dir, "/Trees/", x, ".newick")), simplify = F, USE.NAMES = T)

settings = c("non-sequential", "sequential")
paths = c("non-sequential" = paste0(dir, "/TiDe/baseline/inferenceOutput/"),
          "sequential" = paste0(dir, "/Typewriter/sequential_betterPrior/inferenceOutput/"))
mccTrees = sapply(paths, function(p) {
  sapply(tree_ids, function(x) read.beast(paste0(p, x, ".mcc.tree")), simplify = F, USE.NAMES = T)
}, simplify = F, USE.NAMES = T)
#posTrees = sapply(tree_ids, function(x) read.beast(paste0(path, x, ".thinned.trees")), simplify = F, USE.NAMES = T)

# origin/ experiment duration
Texp = 40

plot_trees <- function(x) {
  
  true = trueTrees[[x]]
  mcc_ns = mccTrees[["non-sequential"]][[x]]
  mcc_s = mccTrees[["sequential"]][[x]]
  
  g_true = ggtree(true, root.position = Texp - max(node.depth.edgelength(true))) +
    geom_tiplab(size = 1) +
    expand_limits(x = c(0, 40))
  
  mcc_ns@phylo = rotateConstr(mcc_ns@phylo, rev(get_taxa_name(g_true)))
  mcc_s@phylo = rotateConstr(mcc_s@phylo, rev(get_taxa_name(g_true)))
  
  g = lapply(list(mcc_ns, mcc_s), function(mcc) {
    ggtree(mcc, root.position = Texp - max(node.depth.edgelength(mcc@phylo))) +
      geom_range(range = 'height_0.95_HPD', center = 'height_median', color = 'lightblue', alpha = 0.5, size = 1.5) +
      geom_nodepoint(aes(color = posterior), size = 2) +
      geom_tiplab(size = 1) +
      expand_limits(x = c(0, 40)) +
      scale_color_continuous(name = "posterior\nsupport", limits = c(0,1), breaks = c(0, 0.5, 1), low = "darkred", high = "darkgreen") +
      theme_tree2() +
      theme(legend.position = c(0.1, 0.8), legend.text = element_text(size = 6), legend.title = element_text(size = 8), 
            legend.key.size = unit(0.4, "cm"), legend.spacing.y = unit(0.1, "cm"))
  })
  
  g = g_true / g[1] / g[2] &
    theme(text = element_text(family = 'lmroman'))
  return(g)
}

# combined plot
g = lapply(tree_ids, plot_trees)
wrap_plots(g, ncol = length(tree_ids), guides = "collect")

ggsave("pdf/SuppFigure_seqTrees.pdf", height = 10, width = 15)
