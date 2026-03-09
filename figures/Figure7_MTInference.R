# Multi-type inference outcomes

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(showtext)
library(patchwork)
library(HDInterval)

# plotting options
palette_types = c("0" = "#56B4E9", "1" = "#009E73", "2" = "#E69F00", "3" = "#CC79A7")
palette_trans = c("distinct" = "#56B4E9", "hierarchical" = "#D55E00")
labels_trans = c("distinct" = "terminal", "hierarchical" = "chain-like")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman') +
            theme(plot.title.position = "plot", plot.subtitle = element_text(margin = margin(l = 30))))

# paths
data_dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/heterog"
code_dir = "~/Projects/celldev"
setwd(file.path(code_dir, "figures"))

# objects
transitions = c("distinct", "hierarchical")

# seeds to keep for results
seeds = sapply(transitions, function(t) {
  trees = read.csv(paste0(code_dir, "/heterog/simulation_trees/tree_params_", t, ".csv")) 
  return(trees$seed) 
}) %>% as.data.frame

nSim = 20

# prior for migration rate
prior = c(hdi(qlnorm, 0.95, meanlog = -3.25, sdlog = 1),
          "median" = qlnorm(0.5, meanlog = -3.25, sdlog = 1))


plot_migration <- function(model, transition, title) {

  infO = readRDS(file.path(data_dir, model, transition, 'analysisOutput', 'infOutput.Rdat')) 

  data = infO %>%
    filter(seed %in% seeds[[transition]]) %>%
    mutate(ID = c(1:nrow(.))) %>%
    select(ID, starts_with('migrationRate')) %>%
    pivot_longer(cols = starts_with("migrationRate"),
                 names_to = c("type", "stat"),
                 names_pattern = "migrationRate\\.t(\\d+)_(\\w+)") %>%
    pivot_wider(names_from = stat, values_from = value)
  
  if (transition == "distinct") {
    data = data %>%
      mutate(true = case_when(type == "1" ~ 0.04,
                              type == "2" ~ 0.03,
                              type == "3" ~ 0.02))
  } else if (transition == "hierarchical") {
    data = data %>%
      mutate(true = case_when(type == "1" ~ 0.08,
                              type == "2" ~ 0.06,
                              type == "3" ~ 0.04))
  }
  
  g = ggplot(data, aes(x = ID, y = median, color = type)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.1, alpha = 0.6, linewidth = 1, show.legend = F) + 
    geom_hline(aes(yintercept = true), color = 'black', alpha = 0.4) +
    geom_hline(aes(yintercept =  prior["lower"]), color = 'black', alpha = 0.4, linetype = 'dotted') +
    geom_hline(aes(yintercept =  prior["upper"]), color = 'black', alpha = 0.4, linetype = 'dotted') +
    scale_y_continuous(limits = c(0,1), breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette_types) +
    facet_wrap(vars(type), labeller = as_labeller(function(x) paste("to type", x))) +
    labs(subtitle = title, x = "Simulation", y = "Posterior transition rate") 
  
  return(g)
}


plot_mcc <- function(model) {
  mccO = sapply(transitions, function(x) {
    data = readRDS(file.path(data_dir, model, x, 'analysisOutput', 'mccOutput.Rdat'))
    data = data %>% filter(seed %in% seeds[[x]])
    }, simplify = F, USE.NAMES = T) %>%
    bind_rows(.id = 'transition')
  
  data = mccO %>%
    select(transition, seed, wRF, prop_type_correct) %>%
    mutate(prop_type_correct = prop_type_correct * 100)
  
  g1 = ggplot(data, aes(x = transition, y = wRF, color = transition)) + 
    geom_boxplot() + #show.legend = FALSE
    geom_hline(yintercept = 0.2, linetype = 'dashed') + # threshold at 0.2
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette_trans, labels = c('terminal', 'chain-like')) +
    theme(axis.ticks.x = element_blank()) + # remove x axis ticks
    labs(x = NULL, y = "Weighted RF distance", color = 'Transitions')
  
  g2 = ggplot(data, aes(x = transition, y = prop_type_correct, color = transition)) + 
    geom_boxplot() + #show.legend = FALSE
    geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 0.2
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_manual(values = palette_trans, labels = labels_trans) +
    theme(axis.ticks.x = element_blank()) + # remove x axis ticks
    labs(x = NULL, y = "Correct ancestral cell types (%)", color = 'Transitions') 
  
  return(g1 + g2)
}
  

p1 = plot_migration('TiDe', 'distinct', labels_trans['distinct'])
p2 = plot_migration('TiDe', 'hierarchical', labels_trans['hierarchical'])
p3 = plot_migration('Typewriter', 'distinct', labels_trans['distinct'])
p4 = plot_migration('Typewriter', 'hierarchical', labels_trans['hierarchical'])
p5 = plot_mcc('TiDe')
p6 = plot_mcc('Typewriter')

p_ns = ((p1 + ggtitle("A: non-sequential recordings")) / p2 / p5) + plot_layout(axis_titles = "collect")
p_s = (p3 + ggtitle("B: sequential recordings")) / p4 / p6 + plot_layout(axis_titles = "collect")

(p_ns | p_s) + plot_layout(guides = "collect", axis_titles = "collect") &
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
          legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 8))
ggsave('pdf/Figure7_MTInference.pdf', height = 10, width = 9)

