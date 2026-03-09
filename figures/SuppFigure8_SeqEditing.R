# Comparison of non- vs. sequential editing

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

# data
paths = c('non-sequential' = file.path(dir, 'TiDe/baseline/analysisOutput'),
          'sequential' = file.path(dir, 'Typewriter/sequential_betterPrior/analysisOutput'))

simStats = sapply(paths, function(p) readRDS(file.path(p, 'simStats.Rdat')), simplify = F, USE.NAMES = T)
infP = sapply(paths, function(p) readRDS(file.path(p, 'infPerformance.Rdat')), simplify = F, USE.NAMES = T)
mccP = sapply(paths, function(p) readRDS(file.path(p, 'mccOutput.Rdat')), simplify = F, USE.NAMES = T)
hpdTrees = sapply(paths, function(p) readRDS(file.path(p, 'hpdTreesOutput.Rdat')), simplify = F, USE.NAMES = T)

# optional: simulation statistics
plot_diversity <- function(simStats) {
  data = lapply(simStats, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, barcodeDiv, saturation)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = barcodeDiv, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Prop. unique barcodes (%)")
  
  return(g)
}


plot_bias <- function(infP) {
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = bias, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.1) + # line at 0
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias") +
    facet_grid(rows = vars(parameter), labeller = labeller(parameter = params_labels), scales = "free_y")
  
  return(g)
}


plot_hpd_width <- function(infP) {
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = hpd_width, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_x_discrete(labels = params_labels) +
    expand_limits(y = c(0, 1)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative HPD width", color = NULL) +
    facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
  
  return(g)
}


plot_hpd_prop <- function(infP) {
  params = c('birthRate', 'deathRate')
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_proportion') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_prop') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
  g = ggplot(data, aes(x = setting, y = hpd_prop, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_x_discrete(labels = params_labels) +
    expand_limits(y = 0) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "HPD proportion", color = NULL) +
    facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
  
  return(g)
}


plot_wRF <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, wRF)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = wRF, color = tree)) + 
    geom_boxplot(outlier.size = 0.5) + # with legend
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Weighted RF distance", color = NULL)
  
  return(g)
}


plot_shPI <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, shPI)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = shPI, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared PI", color = NULL)
  
  return(g)
}


plot_KS <- function(mccP) {
  data = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, KS)}) %>%
    bind_rows(.id = 'setting') %>%
    mutate(tree = factor(tree, levels = trees))
  
  g = ggplot(data, aes(x = setting, y = KS, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "KS distance", color = NULL)
  
  return(g)
}

p1 = plot_wRF(mccP)
p2 = plot_shPI(mccP)
p3 = plot_KS(mccP)
p4 = plot_bias(infP)
p5 = plot_hpd_width(infP)


(p1 + ggtitle('A') + p2 + p3) /
  (p4 + ggtitle('B') + p5) +
  plot_layout(guides = 'collect', heights = c(1, 2.6)) &
  theme(plot.title.position = "plot",
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        legend.position = "bottom", legend.direction = "vertical", legend.justification = "right")
ggsave('pdf/SuppFigure8_SeqEditing.pdf', height = 12, width = 9)

