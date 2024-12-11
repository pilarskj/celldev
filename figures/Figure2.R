# Figure 2

library(dplyr)
library(tidyr)
library(ggplot2)
library(showtext)
library(patchwork)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman'))


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


plot_coverage_lollipop <- function(infP) {
  data = lapply(infP, '[[', 'coverage') %>% 
    bind_rows(.id = 'tree') %>% 
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'coverage') %>% 
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = params)) %>%
    group_by(tree, parameter) %>% 
    summarise(prop = sum(coverage) * 100 / nSim)
  
  g = ggplot(data, aes(x = parameter, y = prop, color = tree)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.2) +
    geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 80%
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Coverage (%)", color = NULL) 

  return(g)
} # add tree topology?


plot_bias_violin <- function(infP) {
  data = lapply(infP, '[[', 'bias_rel') %>% 
    bind_rows(.id = 'tree') %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = params))
  
  g = ggplot(data, aes(x = parameter, y = bias, color = tree)) +
    geom_violin(position = position_dodge(0.6), show.legend = FALSE) +
    # stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 0.5), geom = "pointrange", size = 0.1, position = position_dodge(0.6)) +
    geom_hline(yintercept = 0, linetype = 'dashed') + # line at 0
    ylim(-1, 2) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias", color = NULL)

  return(g)
}


plot_hpd_violin <- function(infP) {
  par = params[1:3]
  par_labels = params_labels[1:3]
  data = lapply(infP, '[[', 'hpd_proportion') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, all_of(par)) %>%
    pivot_longer(cols = all_of(par), names_to = 'parameter', values_to = 'hpd_proportion') %>%
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = par))

  g = ggplot(data, aes(x = parameter, y = hpd_proportion, color = tree)) +
    geom_violin(position = position_dodge(0.6), show.legend = FALSE) +
    geom_hline(yintercept = 0.1, linetype = 'dashed') + # threshold at 0.2
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    scale_x_discrete(labels = par_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "HPD proportion", color = NULL)
  
  return(g)
}


plot_wRF_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, wRF)
  
  g = ggplot(data, aes(x = tree, y = wRF, color = tree)) + # x = metric, y = value
    geom_boxplot(show.legend = FALSE) + # position = position_dodge(0.6)
    geom_hline(yintercept = 0.2, linetype = 'dashed') + # threshold at 0.2
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    theme(axis.ticks.x = element_blank()) + # remove x axis ticks
    labs(x = NULL, y = "Weighted RF distance", color = NULL)

  return(g)
}


plot_summary <- function(path_to_res, title) {
  infP = readRDS(file.path(path_to_res, "infPerformance.Rdat"))
  mccP = readRDS(file.path(path_to_res, "mccOutput.Rdat"))
  
  g1 = plot_coverage_lollipop(infP)
  g2 = plot_bias_violin(infP)
  g3 = plot_hpd_violin(infP)
  g4 = plot_wRF_boxplot(mccP)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 + g4 +
    plot_layout(ncol = 4, nrow = 1, widths = c(2, 2, 1.3, 1))
  return(g)
}

g_td = plot_summary(path_td, 'A: non-sequential recordings') 
g_tw = plot_summary(path_tw, 'B: sequential recordings')

png('Figure2.png', height = 8, width = 11, units = "in", res = 300)
g_td / g_tw +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))
dev.off()

