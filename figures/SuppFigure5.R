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

# optional
plot_diversity_boxplot <- function(simStats) {
  data = lapply(simStats, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, barcodeDiv)}) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = barcodeDiv)) + 
    geom_boxplot(outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Barcode diversity")
  
  return(g)
}


plot_bias_violin <- function(infP) {
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
      mutate(parameter = factor(parameter, levels = params)) }) %>%
    bind_rows(.id = 'setting') 
  
  # g = ggplot(data, aes(x = setting, y = bias)) + 
  #   geom_violin() + 
  #   labs(x = NULL, y = "Relative bias") +
  #   facet_grid(rows = vars(parameter), labeller = as_labeller(params_labels))
  
  g = ggplot(data, aes(x = parameter, y = bias, color = setting)) + 
    geom_violin(show.legend = FALSE) + 
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette[c(2,5)]) +
    labs(x = NULL, y = "Relative bias") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  
  return(g)
}


plot_hpd_violin <- function(infP) {
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params)) }) %>%
    bind_rows(.id = 'setting') 
  
  # g = ggplot(data, aes(x = setting, y = hpd_width)) + 
  #   geom_violin() + 
  #   expand_limits(y = c(0, 1)) +
  #   scale_y_continuous(breaks = scales::pretty_breaks()) +
  #   labs(x = NULL, y = "Relative HPD width") +
  #   facet_grid(rows = vars(parameter), labeller = as_labeller(params_labels), scales = 'free_y')
  
  g = ggplot(data, aes(x = parameter, y = hpd_width, color = setting)) + 
    geom_violin() + 
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette[c(2,5)]) +
    labs(x = NULL, y = "Relative HPD width", color = NULL) +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.9),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))
  
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


plot_Nye_boxplot <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, Nye)}) %>%
    bind_rows(.id = 'setting') 

  g = ggplot(data, aes(x = setting, y = Nye)) + 
    geom_boxplot(outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Nye similarity")
  
  return(g)
}


plot_shPI_boxplot <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, shPI)}) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = shPI)) + 
    geom_boxplot(outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Shared Phylogenetic Information")
  
  return(g)
}


p1 = plot_bias_violin(infP)
p2 = plot_hpd_violin(infP)
p3 = plot_Nye_boxplot(mccP)
p4 = plot_shPI_boxplot(mccP)

png('SuppFigure5.png', height = 8, width = 10, units = "in", res = 300)
(p1 + ggtitle('A') + p2) /
   (p3 + ggtitle('B') + p4) #+
#   plot_layout(heights = c(3, 1))
dev.off()
