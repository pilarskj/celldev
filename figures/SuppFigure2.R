# Supplementary Figure 2

library(dplyr)
library(tidyr)
library(ggplot2)
library(showtext)
library(patchwork)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman') +
            theme(axis.ticks.x = element_blank()))


# paths to infPerformance (TiDe and Typewriter baseline)
path_td = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1/TiDe/baseline/analysisOutput"
path_tw = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1/Typewriter/baseline/analysisOutput"

# parameters and labels
nSim = 20
params = c('birthRate', 'deathRate', 'editRate', 'treeHeight', 'treeLength')
params_labels = c('birthRate' = 'birth rate', 
                  'deathRate' = 'death rate', 
                  'editRate' = 'editing rate', 
                  'treeHeight' = 'tree height', 
                  'treeLength' = 'tree length')

trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")
tree_labels = c("tree_s" = "synchronous trees",
                "tree_ss" = "synchronous trees with sampling",
                "tree_sd" = "synchronous trees with cell death",
                "tree_sds" = "synchronous trees with cell death and sampling",
                "tree_bd" = "birth-death trees with sampling")


# plots for HPD trees
plot_coverage_lollipop <- function(hpdTrees) {
  data = hpdTrees %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    group_by(tree) %>% 
    summarise(prop = sum(recovery) * 100 / nSim)
  
  g = ggplot(data, aes(x = tree, y = prop, color = tree)) +
    geom_point(size = 2) +
    geom_linerange(aes(ymin = 0, ymax = prop), linewidth = 0.2) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = NULL) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Tree coverage (%)", color = NULL) 
  
  return(g)
}


plot_size_boxplot <- function(hpdTrees) {
  data = hpdTrees %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees),
           size = nTrees / nTotal) %>%
    select(tree, seed, size)
  
  g = ggplot(data, aes(x = tree, y = size, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_discrete(labels = NULL) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative size of\n95% credible tree set", color = NULL) 
  
  return(g)
}


# alternative metrics for MCC trees
plot_Nye_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, Nye)
  
  g = ggplot(data, aes(x = tree, y = Nye, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) + 
    geom_hline(yintercept = 0.8, linetype = 'dashed') + # threshold at 0.8
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Nye similarity", color = NULL)
  
  return(g)
}


plot_shPI_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, shPI)
  
  g = ggplot(data, aes(x = tree, y = shPI, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) + 
    geom_hline(yintercept = 0.8, linetype = 'dashed') + # threshold at 0.8
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared Phylogenetic\nInformation", color = NULL)
  
  return(g)
}



plot_summary <- function(path_to_res, title) {
  hpdTrees = readRDS(file.path(path_to_res, "hpdTreesOutput.Rdat"))
  mccP = readRDS(file.path(path_to_res, "mccOutput.Rdat"))
  
  g1 = plot_coverage_lollipop(hpdTrees)
  g2 = plot_size_boxplot(hpdTrees)
  g3 = plot_Nye_boxplot(mccP)
  g4 = plot_shPI_boxplot(mccP)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 + g4 +
    plot_layout(ncol = 4, nrow = 1)
  return(g)
}

g_td = plot_summary(path_td, 'A: non-sequential recordings') 
g_tw = plot_summary(path_tw, 'B: sequential recordings')

pdf('figures/SuppFigure2.pdf', height = 8, width = 11)
g_td / g_tw +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")
dev.off()

