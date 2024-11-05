# Figure 3

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
dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/Part1"
paths = c('non-sequential' = file.path(dir, 'TiDe/baseline/analysisOutput'),
          'sequential' = file.path(dir, 'Typewriter/sequential_betterPrior/analysisOutput'))

mccP = sapply(paths, function(p) readRDS(file.path(p, 'mccOutput.Rdat')), simplify = F, USE.NAMES = T)
infP = sapply(paths, function(p) readRDS(file.path(p, 'infPerformance.Rdat')), simplify = F, USE.NAMES = T)

plot_wRF_boxplot <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, wRF)}) %>%
    bind_rows(.id = 'setting') 

  g = ggplot(data, aes(x = setting, y = wRF)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 't.test', map_signif_level = T) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Weighted RF distance")
    
  return(g)
}


plot_hpd_boxplot <- function(infP) {
  params = c('birthRate', 'deathRate')
  params_labels = c('birthRate' = 'birth rate', 'deathRate' = 'death rate')
  
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_proportion') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_proportion') %>%
      mutate(parameter = factor(parameter, levels = params)) }) %>%
    bind_rows(.id = 'setting') 

  g = ggplot(data, aes(x = setting, y = hpd_proportion)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 't.test', map_signif_level = T) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    labs(x = NULL, y = "HPD proportion") +
    facet_wrap(vars(parameter), labeller = as_labeller(params_labels))
  
  return(g)
}


png('Figure4.png', height = 5, width = 8, units = "in", res = 300)
plot_wRF_boxplot(mccP) + plot_hpd_boxplot(infP) + plot_layout(widths = c(1, 2)) + plot_annotation(tag_levels = 'A')
dev.off()


