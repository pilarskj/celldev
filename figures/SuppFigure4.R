# Supplementary Figure 4

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(showtext)
library(patchwork)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman'))

# number of simulations per setting
nSim = 20

# tree labels
trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")
tree_labels = c("tree_s" = "synchronous trees",
                "tree_ss" = "synchronous trees with sampling",
                "tree_sd" = "synchronous trees with cell death",
                "tree_sds" = "synchronous trees with cell death and sampling",
                "tree_bd" = "birth-death trees with sampling")

# settings
settings_td = data.frame(
  setting = c('baseline', 'editRate_0.01', 'editRate_0.15', 'nTargets_5', 'nTargets_40', 'scarringDuration_20', 'scarringWindow', 'scarringHeight_20'),
  level = c('', '0.01', '0.15', '5', '40', 'first half', 'mid', 'second half')) %>%
  mutate(manipulation = factor(case_when(str_starts(setting, "editRate") ~ "editing rate",
                                         str_starts(setting, "nTargets") ~ "no. targets",
                                         str_starts(setting, "scarring") ~ "editing window",
                                         .default = "baseline"),
                               levels = c("baseline", "editing rate", "no. targets", "editing window")))
  
settings_tw = data.frame(
  setting = c('baseline', 'editRate_0.01', 'editRate_0.15', 'nTapes_5', 'nTapes_40', 'tapeLength_2', 'tapeLength_10'),
  level = c('', '0.01', '0.15', '5', '40', '2', '10')) %>%
  mutate(manipulation = factor(case_when(str_starts(setting, "editRate") ~ "editing rate",
                                         str_starts(setting, "nTapes") ~ "no. tapes",
                                         str_starts(setting, "tapeLength") ~ "tape length",
                                         .default = "baseline"),
                               levels = c("baseline", "editing rate", "no. tapes", "tape length")))

# paths
path_td = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1/TiDe"
path_tw = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1/Typewriter"


plot_treecov_lollipop <- function(hpdTrees, settings) {
  
  data = sapply(hpdTrees, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      group_by(tree) %>%
      summarise(prop = sum(recovery) * 100 / nSim) %>%
      select(tree, prop)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(prop) %>% median
  
  g = ggplot(data, aes(x = level, y = prop, color = tree)) + 
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.2) +
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Tree coverage (%)", color = NULL) +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_treesize_boxplot <- function(hpdTrees, settings) {
  
  data = sapply(hpdTrees, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      mutate(size = nTrees / nTotal)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(size) %>% median
  
  g = ggplot(data, aes(x = level, y = size, color = tree)) +
    geom_boxplot(outlier.size = 0.5, show.legend = F) +
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative size of\n95% credible tree set") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_wRF_boxplot <- function(mccP, settings) {
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, wRF)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(wRF) %>% median
  
  g = ggplot(data, aes(x = level, y = wRF, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Weighted RF distance") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_Nye_boxplot <- function(mccP, settings) {
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, Nye)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(Nye) %>% median
  
  g = ggplot(data, aes(x = level, y = Nye, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Nye similarity") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_shPI_boxplot <- function(mccP, settings) {
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, shPI)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(shPI) %>% median
  
  g = ggplot(data, aes(x = level, y = shPI, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared Phylogenetic\nInformation") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_summary <- function(path, settings, title) {
  
  hpdTrees = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'hpdTreesOutput.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  mccP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'mccOutput.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  g1 = plot_wRF_boxplot(mccP, settings)
  g2 = plot_Nye_boxplot(mccP, settings)
  g3 = plot_shPI_boxplot(mccP, settings)
  g4 = plot_treecov_lollipop(hpdTrees, settings)
  g5 = plot_treesize_boxplot(hpdTrees, settings)
  
  
  g = ((g1 + ggtitle(title) + theme(plot.title.position = "plot")) / g2 / g3) | 
    (g4 + g5 + guide_area() + plot_layout(ncol = 1, guides = 'collect'))
  
  return(g)
}


p_td = plot_summary(path_td, settings_td, 'A: non-sequential recordings')
p_tw = plot_summary(path_tw, settings_tw, 'B: sequential recordings')

pdf('figures/SuppFigure4.pdf', height = 12, width = 12) 
p_td
p_tw
dev.off()