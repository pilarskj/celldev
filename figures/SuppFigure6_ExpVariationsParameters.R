# Parameter inference at varying experimental parameters

setwd("~/Projects/celldev/figures")
source("homog_functions.R")
theme_update(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))

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

path_td = file.path(dir, "TiDe")
path_tw = file.path(dir, "Typewriter")


plot_bias <- function(infP, settings) {
  
  data = sapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
      mutate(parameter = factor(parameter, levels = params))}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% group_by(parameter) %>% summarize(median = median(bias))
  
  g = ggplot(data, aes(x = level, y = bias, color = tree)) + 
    geom_boxplot(outlier.size = 0.5) + 
    geom_hline(data = baseline, aes(yintercept = median), linetype = 'dashed') +
    expand_limits(y = c(-0.5, 0.5)) +
    scale_y_continuous(breaks = seq(-1, 2, 0.5)) + 
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias", color = NULL) +
    facet_grid(rows = vars(parameter), cols = vars(manipulation), labeller = labeller(parameter = params_labels),
               scales = 'free', space = 'free_x')
  
  return(g)
}


plot_hpd <- function(infP, settings) {
  
  data = sapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params))}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees))  
  
  baseline = data %>% filter(setting == 'baseline') %>% group_by(parameter) %>% summarize(median = median(hpd_width))
  
  g = ggplot(data, aes(x = level, y = hpd_width, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(data = baseline, aes(yintercept = median), linetype = 'dashed') +
    expand_limits(y = c(0, 1)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    scale_color_manual(values = palette) +
    labs(x = NULL, y = "Relative HPD width") +
    facet_grid(rows = vars(parameter), cols = vars(manipulation), labeller = labeller(parameter = params_labels),
               scales = 'free', space = 'free_x') 
  
  return(g)
}


plot_summary <- function(path, settings, title) {
  
  infP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'infPerformance.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  g1 = plot_bias(infP, settings)
  g2 = plot_hpd(infP, settings)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") + g2 +
    plot_layout(ncol = 1)
  
  return(g)
}


p_td = plot_summary(path_td, settings_td, 'A: non-sequential recordings')
p_tw = plot_summary(path_tw, settings_tw, 'B: sequential recordings')

pdf('SuppFigure6_ExpVariationsParameters.pdf', height = 12, width = 10) 
p_td
p_tw
dev.off()
