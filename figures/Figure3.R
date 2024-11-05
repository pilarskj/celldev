# Figure 3

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
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


plot_wRF_boxplot <- function(path, settings) {
  mccP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'mccOutput.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, wRF)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(wRF) %>% median

  g = ggplot(data, aes(x = level, y = wRF)) + 
    geom_boxplot() + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Weighted RF distance") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
    
  return(g)
}


plot_hpd_pointrage <- function(path, settings) {
  
  infP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'infPerformance.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  data = sapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params))}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level)) 
  
  # just plot editing rate
  data = data %>%
    filter(manipulation %in% c('baseline', 'editing rate')) %>%
    mutate(level = ifelse(manipulation == 'baseline', '0.05', as.character(level)))
  
  g = ggplot(data, aes(x = parameter, y = hpd_width, color = level)) +
    stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), geom = "pointrange", size = 0.1, position = position_dodge(0.5)) +
    ylim(0, 4) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette[1:3]) +
    labs(x = NULL, y = "Relative HPD width", color = 'editing rate') +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.9),
          legend.background = element_rect(linetype = "dashed", color = 'black', linewidth = 0.2),
          legend.title = element_text(size = 10), legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), axis.text.y = element_text(size = 8))
  
  return(g)
}


pA = plot_wRF_boxplot(path_td, settings_td)
pB = plot_hpd_pointrage(path_td, settings_td)
pC = plot_wRF_boxplot(path_tw, settings_tw)
pD = plot_hpd_pointrage(path_tw, settings_tw)

png('Figure3.png', height = 8, width = 9, units = "in", res = 300) 
p_td = pA + ggtitle('A: non-sequential recordings') + theme(plot.title.position = "plot") + pB + plot_layout(ncol = 2, widths = c(2, 1))
p_tw = pC + ggtitle('B: sequential recordings') + theme(plot.title.position = "plot") + pD + plot_layout(ncol = 2, widths = c(2, 1))
p_td / p_tw  
dev.off()


