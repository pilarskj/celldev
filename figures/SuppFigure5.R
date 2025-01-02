# Supplementary Figure 5

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(scales)
library(showtext)
library(patchwork)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman') + theme(plot.title.position = "plot"))

# data
dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1"
paths = c('non-sequential' = file.path(dir, 'TiDe/baseline/analysisOutput'),
          'sequential' = file.path(dir, 'Typewriter/sequential_betterPrior/analysisOutput'))

simStats = sapply(paths, function(p) readRDS(file.path(p, 'simStats.Rdat')), simplify = F, USE.NAMES = T)
infP = sapply(paths, function(p) readRDS(file.path(p, 'infPerformance.Rdat')), simplify = F, USE.NAMES = T)
mccP = sapply(paths, function(p) readRDS(file.path(p, 'mccOutput.Rdat')), simplify = F, USE.NAMES = T)
hpdTrees = sapply(paths, function(p) readRDS(file.path(p, 'hpdTreesOutput.Rdat')), simplify = F, USE.NAMES = T)

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

# optional
plot_diversity_boxplot <- function(simStats) {
  data = lapply(simStats, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, barcodeDiv)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = barcodeDiv, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Barcode diversity")
  
  return(g)
}


plot_bias_violin <- function(infP) {
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = bias, color = tree)) + 
    geom_violin(show.legend = F, position = position_dodge(0.6)) + 
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias") +
    facet_grid(rows = vars(parameter), labeller = labeller(parameter = params_labels))
  
  return(g)
}


plot_hpd_violin <- function(infP) {
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = hpd_width, color = tree)) + 
    geom_violin(position = position_dodge(0.6)) + 
    scale_x_discrete(labels = params_labels) +
    expand_limits(y = 1) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative HPD width", color = NULL) +
    facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
  
  return(g)
}

# boring
plot_treecov_lollipop <- function(hpdTrees) {
  data = lapply(hpdTrees, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      summarise(prop = sum(recovery) * 100 / 100) %>%
      select(prop)}) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = prop)) + 
    geom_point(size = 2) +
    geom_linerange(aes(ymin = 0, ymax = prop), linewidth = 0.2) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(x = NULL, y = "Tree coverage (%)", color = NULL)
  
  return(g)
}

# boring
plot_treesize_boxplot <- function(hpdTrees) {
  data = lapply(hpdTrees, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      mutate(size = nTrees / nTotal) %>%
      select(tree, seed, size)}) %>%
    bind_rows(.id = 'setting')
  
  g = ggplot(data, aes(x = setting, y = size)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Relative size of\n95% credible tree set")
  
  return(g)
}


plot_wRF_boxplot <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, wRF)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = wRF, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Weighted RF distance")
  
  return(g)
}


plot_Nye_boxplot <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, Nye)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))

  g = ggplot(data, aes(x = setting, y = Nye, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Nye similarity")
  
  return(g)
}


plot_shPI_boxplot <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, shPI)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = shPI, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared Phylogenetic Information")
  
  return(g)
}


p1 = plot_bias_violin(infP)
p2 = plot_hpd_violin(infP)
p3 = plot_wRF_boxplot(mccP)
p4 = plot_Nye_boxplot(mccP)
p5 = plot_shPI_boxplot(mccP)


pdf('figures/SuppFigure5.pdf', height = 12, width = 10)
(p1 + ggtitle('A') + p2) /
  (p3 + ggtitle('B') + p4 + p5) +
  plot_layout(guides = 'collect', heights = c(2, 1)) & 
  theme(legend.position = "bottom", legend.direction = "vertical", legend.justification = "right")
dev.off()
