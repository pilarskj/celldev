# Figure 5

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggsignif)
library(showtext)
library(patchwork)
library(HDInterval)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman'))

# data
dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part2"

# objects
transitions = c("distinct", "hierarchical")

# seeds to keep for results
seeds = sapply(transitions, function(t) {
  trees = read.csv(paste0("~/Projects/celldev/heterog/simulation_trees/tree_params_", t, ".csv")) 
  return(trees$seed) 
}) %>% as.data.frame

nSim = 20

# prior for migration rate
prior = c(hdi(qlnorm, 0.95, meanlog = -3.25, sdlog = 1),
          "median" = qlnorm(0.5, meanlog = -3.25, sdlog = 1))


plot_migration <- function(model, transition, title) {

  infO = readRDS(file.path(dir, model, transition, 'analysisOutput', 'infOutput.Rdat')) 

  data = infO %>%
    filter(seed %in% seeds[[transition]]) %>%
    mutate(ID = c(1:nrow(.))) %>%
    select(ID, starts_with('migrationRate')) %>%
    pivot_longer(cols = starts_with("migrationRate"),
                 names_to = c("t", "stat"),
                 names_pattern = "migrationRate\\.t(\\d+)_(\\w+)") %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(type = paste0("to type ", t))
  
  if (transition == "distinct") {
    data = data %>%
      mutate(true = case_when(t == "1" ~ 0.04,
                              t == "2" ~ 0.03,
                              t == "3" ~ 0.02))
  } else if (transition == "hierarchical") {
    data = data %>%
      mutate(true = case_when(t == "1" ~ 0.08,
                              t == "2" ~ 0.06,
                              t == "3" ~ 0.04))
  }
  
  g = ggplot(data, aes(x = ID, y = median, color = type)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.1, alpha = 0.6, linewidth = 1, show.legend = F) + 
    geom_hline(aes(yintercept = true), color = 'black', alpha = 0.4) +
    geom_hline(aes(yintercept =  prior["lower"]), color = 'black', alpha = 0.4, linetype = 'dotted') +
    geom_hline(aes(yintercept =  prior["upper"]), color = 'black', alpha = 0.4, linetype = 'dotted') +
    ylim(c(0, 1)) + 
    scale_color_manual(values = palette[3:5]) +
    facet_wrap(vars(type)) +
    labs(subtitle = title, x = "Simulation", y = "Posterior migration rate") 
  
  return(g)
}


plot_mcc <- function(model) {
  mccO = sapply(transitions, function(x) {
    data = readRDS(file.path(dir, model, x, 'analysisOutput', 'mccOutput.Rdat'))
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
    scale_color_manual(values = palette[1:2], labels = c('terminal', 'chain-like')) +
    theme(axis.ticks.x = element_blank()) + # remove x axis ticks
    labs(x = NULL, y = "Weighted RF distance", color = 'transitions')
  
  g2 = ggplot(data, aes(x = transition, y = prop_type_correct, color = transition)) + 
    geom_boxplot() + #show.legend = FALSE
    geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 0.2
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_manual(values = palette[1:2], labels = c('terminal', 'chain-like')) +
    theme(axis.ticks.x = element_blank()) + # remove x axis ticks
    labs(x = NULL, y = "Correct types (%)", color = 'transitions') 
  
  return(g1 + g2)
}
  

p1 = plot_migration('TiDe', 'distinct', 'terminal transitions')
p2 = plot_migration('TiDe', 'hierarchical', 'chain-like transitions')
p3 = plot_migration('Typewriter', 'distinct', 'terminal transitions')
p4 = plot_migration('Typewriter', 'hierarchical', 'chain-like transitions')
pA = p1 + ggtitle('non-sequential recordings') + p2 + 
  p3 + ggtitle('sequential recordings') + p4 #+ 
#plot_layout(axis_titles = 'collect')

p5 = plot_mcc('TiDe')
p6 = plot_mcc('Typewriter')
pB = p5 / p6 + plot_layout(guides = 'collect') &
  theme(legend.position = "bottom",
        legend.background = element_rect(linetype = "dashed", color = 'black', linewidth = 0.2),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))


pdf('Figure5.pdf', height = 8, width = 12)#, units = "in", res = 300) 
(pA | pB) + plot_layout(widths = c(2, 1))
dev.off()




