# Figure 5 Supp # TODO

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(showtext)
library(patchwork)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman'))

# data
dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part2"

# objects
transitions = c("distinct", "hierarchical")
transitions_labels = c("distinct" = "terminal", "hierarchical" = "chain-like")

# seeds to keep for results
seeds = sapply(transitions, function(t) {
  trees = read.csv(paste0("~/Projects/celldev/heterog/simulation_trees/tree_params_", t, ".csv")) 
  return(trees$seed) 
}) %>% as.data.frame

nSim = 20

# inferred parameters
params = c("birthRate.t0", "birthRate.t1", "birthRate.t2", "birthRate.t3",
           "deathRate.t0", "deathRate.t1", "deathRate.t2", "deathRate.t3", 
           "migrationRate.t1", "migrationRate.t2", "migrationRate.t3",
           "editRate", "treeHeight", "treeLength")

params_labels = c("birthRate.t0" = "birth rate (type 0)", "birthRate.t1" = "birth rate (type 1)",
                  "birthRate.t2" = "birth rate (type 2)", "birthRate.t3" = "birth rate (type 3)",
                  "deathRate.t0" = "death rate (type 0)", "deathRate.t1" = "death rate (type 1)", 
                  "deathRate.t2" = "death rate (type 2)", "deathRate.t3" = "death rate (type 3)",
                  "migrationRate.t1" = "migration rate (to type 1)", 
                  "migrationRate.t2" = "migration rate (to type 2)",
                  "migrationRate.t3" = "migration rate (to type 3)",
                  "editRate" = "editing rate", "treeHeight" = "tree height", "treeLength" = "tree length")


plot_coverage_lollipop <- function(paths) {
  data = lapply(infP, '[[', 'coverage') 
  data$distinct = data$distinct %>% filter(seed %in% seeds$distinct)
  data$hierarchical = data$hierarchical %>% filter(seed %in% seeds$hierarchical)
  data = data %>% 
    bind_rows(.id = 'transition') %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'coverage') %>% 
    mutate(parameter = factor(parameter, levels = params)) %>%
    group_by(transition, parameter) %>% 
    summarise(prop = sum(coverage) * 100 / nSim)
  
  g = ggplot(data, aes(x = parameter, y = prop, color = transition)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) +
    geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.2) +
    geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 80%
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = transitions_labels) +
    labs(x = NULL, y = "Coverage (%)")#, color = NULL) 
  
  return(g)
} 


plot_bias_violin <- function(infP) {
  data = lapply(infP, '[[', 'bias_rel') 
  data$distinct = data$distinct %>% filter(seed %in% seeds$distinct)
  data$hierarchical = data$hierarchical %>% filter(seed %in% seeds$hierarchical)
  data = data %>% 
    bind_rows(.id = 'transition') %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
    mutate(parameter = factor(parameter, levels = params))
  
  g = ggplot(data, aes(x = parameter, y = bias, color = transition)) +
    #geom_violin(position = position_dodge(0.6), show.legend = FALSE) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 0.5), geom = "pointrange", size = 0.1, position = position_dodge(0.6)) +
    geom_hline(yintercept = 0, linetype = 'dashed') + # line at 0
    #ylim(-1, 2) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, transitions_labels) +
    labs(x = NULL, y = "Relative bias")#, color = NULL)
  
  return(g)
}



plot_summary <- function(paths, title) {
  infP = lapply(paths, function(x) readRDS(file.path(x, "infPerformance.Rdat")))
  #mccP = readRDS(file.path(path, "mccOutput.Rdat"))
  
  g1 = plot_coverage_lollipop(infP)
  g2 = plot_bias_violin(infP)
  g3 = plot_hpd_violin(infP)
  g4 = plot_wRF_boxplot(mccP)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 + g4 +
    plot_layout(ncol = 4, nrow = 1, widths = c(2, 2, 1.3, 1))
  return(g)
}

paths = sapply(names(transitions), function(x) file.path(dir, 'TiDe', x, 'analysisOutput'))
