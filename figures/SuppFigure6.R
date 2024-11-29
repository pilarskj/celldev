# Supplementary Figure 6

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
theme_set(theme_classic(base_size = 12, base_family = 'lmroman') +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
                  axis.text.y = element_text(size = 8)))

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
inf_params = c("birthRate.t0", "birthRate.t1", "birthRate.t2", "birthRate.t3",
               "deathRate.t0", "deathRate.t1", "deathRate.t2", "deathRate.t3", 
               "migrationRate.t1", "migrationRate.t2", "migrationRate.t3",
               "editRate", "treeHeight", "treeLength")

inf_params_labels = c("birthRate.t0" = "birth rate (type 0)", "birthRate.t1" = "birth rate (type 1)",
                      "birthRate.t2" = "birth rate (type 2)", "birthRate.t3" = "birth rate (type 3)",
                      "deathRate.t0" = "death rate (type 0)", "deathRate.t1" = "death rate (type 1)", 
                      "deathRate.t2" = "death rate (type 2)", "deathRate.t3" = "death rate (type 3)",
                      "migrationRate.t1" = "migration rate (to type 1)", 
                      "migrationRate.t2" = "migration rate (to type 2)",
                      "migrationRate.t3" = "migration rate (to type 3)",
                      "editRate" = "editing rate", "treeHeight" = "tree height", "treeLength" = "tree length")

params = c('birthRate', 'deathRate', 'migrationRate', 'editRate', 'treeHeight', 'treeLength')
params_labels = c('birthRate' = 'birth rate', 
                  'deathRate' = 'death rate', 
                  'migrationRate' = 'migration rate', 
                  'editRate' = 'editing rate', 
                  'treeHeight' = 'tree height', 
                  'treeLength' = 'tree length')
type_labels = c("0" = "type 0", "1" = "(to) type 1", "2" = "(to) type 2", "3" = "(to) type 3")


plot_coverage_lollipop <- function(paths) {
  data = lapply(infP, '[[', 'coverage') 
  data$distinct = data$distinct %>% filter(seed %in% seeds$distinct)
  data$hierarchical = data$hierarchical %>% filter(seed %in% seeds$hierarchical)
  data = data %>%
    bind_rows(.id = 'transition') %>%
    pivot_longer(cols = all_of(inf_params), names_to = 'parameter', values_to = 'coverage') %>%
    group_by(transition, parameter) %>%
    summarise(prop = sum(coverage) * 100 / nSim) %>% 
    mutate(type = str_extract(parameter, '[0-3]'), 
           parameter = factor(str_remove(parameter, '\\.t[0-3]'), levels = params))

  g = ggplot(data, aes(x = parameter, y = prop, color = type)) +
    geom_point(size = 2, position = position_dodge(width = 0.8)) +
    geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.8), linewidth = 0.2) +
    geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 80%
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, limits = c(names(type_labels), ""), 
                       labels = type_labels, na.value = "black") +
    labs(x = NULL, y = "Coverage (%)", color = NULL) +
    facet_wrap(vars(transition), labeller = labeller(transition = transitions_labels))
  return(g)
} 


plot_bias_boxplot <- function(infP) {
  data = lapply(infP, '[[', 'bias_rel') 
  data$distinct = data$distinct %>% filter(seed %in% seeds$distinct)
  data$hierarchical = data$hierarchical %>% filter(seed %in% seeds$hierarchical)
  data = data %>% 
    bind_rows(.id = 'transition') %>%
    pivot_longer(cols = all_of(inf_params), names_to = 'parameter', values_to = 'bias') %>%
    mutate(type = str_extract(parameter, '[0-3]'), 
           parameter = factor(str_remove(parameter, '\\.t[0-3]'), levels = params))
  
  g = ggplot(data, aes(x = parameter, y = bias, color = type)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8), show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = 'dashed') + # line at 0
    scale_x_discrete(labels = params_labels) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, na.value = "black") +
    labs(x = NULL, y = "Relative bias") +
    facet_wrap(vars(transition), labeller = labeller(transition = transitions_labels))
  
  return(g)
}


plot_hpd_boxplot <- function(infP) {
  data = lapply(infP, '[[', 'hpd_proportion') 
  data$distinct = data$distinct %>% filter(seed %in% seeds$distinct)
  data$hierarchical = data$hierarchical %>% filter(seed %in% seeds$hierarchical)
  data = data %>% 
    bind_rows(.id = 'transition') %>%
    pivot_longer(cols = all_of(inf_params), names_to = 'parameter', values_to = 'hpd_proportion') %>%
    mutate(type = str_extract(parameter, '[0-3]'), 
           parameter = factor(str_remove(parameter, '\\.t[0-3]'), levels = params)) %>% 
    filter(parameter %in% c('birthRate', 'deathRate', 'migrationRate'))
  
  g = ggplot(data, aes(x = parameter, y = hpd_proportion, color = type)) +
    geom_boxplot(outlier.size = 0.5, position = position_dodge(0.8, preserve = "single"), show.legend = FALSE) +
    scale_x_discrete(labels = params_labels) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette) +
    labs(x = NULL, y = "HPD proportion") +
    facet_wrap(vars(transition), labeller = labeller(transition = transitions_labels))
  
  return(g)
}


plot_summary <- function(recorder, title) {
  infP = sapply(transitions, function(x) 
    readRDS(file.path(dir, recorder, x, 'analysisOutput', 'infPerformance.Rdat')),
    simplify = F, USE.NAMES = T)
  
  g1 = plot_coverage_lollipop(infP)
  g2 = plot_bias_boxplot(infP)
  g3 = plot_hpd_boxplot(infP)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 +
    plot_layout(ncol = 3, nrow = 1) #, widths = c(2, 2, 1.3, 1))
  return(g)
}

g_td = plot_summary('TiDe', 'A: non-sequential recordings')
g_tw = plot_summary('Typewriter', 'B: sequential recordings')

png('SuppFigure6.png', height = 10, width = 14, units = "in", res = 300) 
g_td / g_tw + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')
dev.off()