# Branch length distributions

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

# all runs
runs = with(expand.grid(tree = trees, sim = seq_len(nSim)), paste0(tree, "_", sim))

plot_branch_dist <- function(trueDir, infDir) {
  
  # load trees (.newick or .mcc.tree format)
  trueTrees = sapply(runs, function(x) read.tree(paste0(trueDir, "/", x, ".newick")), simplify = F, USE.NAMES = T)
  infTrees = sapply(runs, function(x) {
    f = paste0(infDir, "/", x, ".mcc.tree") 
    tree = read.beast(f)
    return(as.phylo(tree))
  }, simplify = F, USE.NAMES = T)
  
  # collect branch lengths
  df = lapply(trueTrees, function(t) data.frame(branch = t$edge.length)) %>%
    bind_rows(.id = "tree") %>% mutate(group = "true")
  df = rbind(df, lapply(infTrees, function(t) data.frame(branch = t$edge.length)) %>% 
               bind_rows(.id = "tree") %>% mutate(group = "inferred"))
  df = df %>%
    mutate(seed = sub("tree_[a-z]+_", "", tree)) %>%
    mutate(tree = sub("_[0-9]+", "", tree)) %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(df, aes(x = branch, y = tree, linetype = group, color = tree)) +
    geom_density_ridges(fill = NA, scale = 1.25) + #jittered_points = TRUE, point_alpha = 0.05, point_size = 0.0001) + 
    scale_x_continuous(limits = c(-2, 40), breaks = seq(0,40,10)) +
    scale_y_discrete(limits = rev(trees)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = "Branch length", y = NULL, linetype = NULL, color = NULL) +
    theme(axis.text.y = element_blank())
  
  return(g)
}

# # Baseline scenarios
#g_ns = plot_branch_dist(trueDir = paste0(dir, "/Trees"), infDir = paste0(dir, "/TiDe/baseline/inferenceOutput"))
#g_s = plot_branch_dist(trueDir = paste0(dir, "/Trees"), infDir = paste0(dir, "/Typewriter/baseline/inferenceOutput"))
# Noisy data
g_ns = plot_branch_dist(trueDir = paste0(dir, "/TiDe/baseline_noise/filteredTrees"), infDir = paste0(dir, "/TiDe/baseline_noise/inferenceOutput"))
g_s = plot_branch_dist(trueDir = paste0(dir, "/Typewriter/baseline_noise/filteredTrees"), infDir = paste0(dir, "/Typewriter/baseline_noise/inferenceOutput"))
(g_ns + ggtitle("A: non-sequential recordings (noisy)")) +
  (g_s + ggtitle("B: sequential recordings (noisy)")) +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(plot.title.position = "plot",
        legend.position = "bottom", legend.direction = "vertical", legend.justification = "right")
#ggsave('pdf/SuppFigure4_BaselineBranches.pdf', height = 8, width = 8)
ggsave('pdf/SuppFigure12_NoiseBranches.pdf', height = 8, width = 8)


# Comparison
# g1 = plot_branch_dist(trueDir = paste0(dir, "/Trees"), infDir = paste0(dir, "/TiDe/nTargets_5/inferenceOutput"))
# g2 = plot_branch_dist(trueDir = paste0(dir, "/TiDe/baseline_noise/filteredTrees"), infDir = paste0(dir, "/TiDe/baseline_noise/inferenceOutput"))
# g3 = plot_branch_dist(trueDir = paste0(dir, "/Trees"), infDir = paste0(dir, "/Typewriter/nTapes_5/inferenceOutput"))
# g4 = plot_branch_dist(trueDir = paste0(dir, "/Typewriter/baseline_noise/filteredTrees"), infDir = paste0(dir, "/Typewriter/baseline_noise/inferenceOutput"))
# ((g1 + ggtitle("A: non-sequential recordings") + labs(subtitle = "5 targets")) +
#   (g2 + labs(subtitle = "filtered targets"))) /
#   ((g3 + ggtitle("B: sequential recordings") + labs(subtitle = "5 tapes")) + 
#   (g4 + labs(subtitle = "filtered tapes"))) +
#   plot_layout(guides = "collect", axis_titles = "collect") & 
#   theme(plot.title.position = "plot", 
#         legend.position = "bottom", legend.direction = "vertical", legend.justification = "right",
#         axis.text.x = element_text(size = 8), axis.title.x = element_text(size = 8))
# ggsave('SuppFigure_Branches_NoiseComp.pdf', height = 12, width = 8)
