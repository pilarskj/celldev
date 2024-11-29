# Supplementary Figure 3

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

# parameters and labels
params = c('birthRate', 'deathRate', 'editRate', 'treeHeight', 'treeLength')
params_labels = c('birthRate' = 'birth rate', 
                  'deathRate' = 'death rate', 
                  'editRate' = 'editing rate', 
                  'treeHeight' = 'tree height', 
                  'treeLength' = 'tree length')

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

path_td = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1/TiDe"
path_tw = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1/Typewriter"


plot_diversity_boxplot <- function(simStats, settings) {
  
  data = sapply(simStats, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, barcodeDiv)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(barcodeDiv) %>% median
  
  g = ggplot(data, aes(x = level, y = barcodeDiv)) + 
    geom_boxplot(outlier.size = 0.5) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Barcode diversity") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_bias_violin <- function(infP, settings) {
  
  data = sapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
      mutate(parameter = factor(parameter, levels = params))}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% group_by(parameter) %>% summarize(median = median(bias))

  g = ggplot(data, aes(x = level, y = bias, color = manipulation)) + 
    geom_violin(show.legend = F) + 
    geom_hline(data = baseline, aes(yintercept = median), linetype = 'dashed') +
    scale_y_continuous(limits = c(-1, 2), breaks = seq(-1, 2, 1)) +
    scale_color_manual(values = palette) +
    labs(x = NULL, y = "Relative bias") +
    facet_grid(rows = vars(parameter), cols = vars(manipulation), labeller = labeller(parameter = params_labels),
               scales = 'free_x', space = 'free_x')
    
  return(g)
}


plot_hpd_violin <- function(infP, settings) {
  
  data = sapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params))}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% group_by(parameter) %>% summarize(median = median(hpd_width))
  
  g = ggplot(data, aes(x = level, y = hpd_width, color = manipulation)) + 
    geom_violin(show.legend = F) + 
    geom_hline(data = baseline, aes(yintercept = median), linetype = 'dashed') +
    expand_limits(y = c(0, 1)) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette) +
    labs(x = NULL, y = "Relative HPD width") +
    facet_grid(rows = vars(parameter), cols = vars(manipulation), labeller = labeller(parameter = params_labels),
               scales = 'free', space = 'free_x') 
  
  return(g)
}


plot_summary <- function(path, settings, title) {
  
  simStats = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'simStats.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  infP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'infPerformance.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  g1 = plot_diversity_boxplot(simStats, settings)
  g2 = plot_bias_violin(infP, settings)
  g3 = plot_hpd_violin(infP, settings)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 +
    plot_layout(ncol = 1, nrow = 3, heights = c(1, 3, 3))
  
  return(g)
}


p_td = plot_summary(path_td, settings_td, 'A: non-sequential recordings')
p_tw = plot_summary(path_tw, settings_tw, 'B: sequential recordings')

png('SuppFigure3.png', height = 18, width = 12, units = "in", res = 300) 
wrap_plots(p_td, p_tw)
dev.off()
