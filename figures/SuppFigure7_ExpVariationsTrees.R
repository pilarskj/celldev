# Tree inference at varying experimental parameters

setwd("~/Projects/celldev/figures")
source("homog_functions.R")
theme_update(axis.ticks.x = element_blank(), axis.text.y = element_text(size = 8))

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


plot_diversity <- function(simStats, settings) {
  
  data = sapply(simStats, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, barcodeDiv)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(propUnique = barcodeDiv * 100) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(propUnique) %>% median
  
  g = ggplot(data, aes(x = level, y = propUnique, color = tree)) + 
    geom_boxplot(outlier.size = 0.5) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Prop. unique barcodes (%)", color = NULL) +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_treecov_lollipop <- function(hpdTrees, settings) {
  
  data = sapply(hpdTrees, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      group_by(tree) %>%
      summarise(prop = sum(recovery) * 100 / nSim) %>%
      select(tree, prop)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(prop) %>% median
  
  g = ggplot(data, aes(x = level, y = prop, color = tree)) + 
    geom_point(size = 2, position = position_dodge(width = 0.6), show.legend = F) +
    geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.2, show.legend = F) +
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Tree coverage (%)", color = NULL) +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_wRF_boxplot <- function(mccP, settings) {
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, wRF)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(wRF) %>% median
  
  g = ggplot(data, aes(x = level, y = wRF, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Weighted RF distance") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_shPI_boxplot <- function(mccP, settings) {
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, shPI)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(shPI) %>% median
  
  g = ggplot(data, aes(x = level, y = shPI, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared PI") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_KS_boxplot <- function(mccP, settings) {
  
  data = sapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, KS)}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) %>%
    mutate(level = factor(level, levels = settings$level),
           tree = factor(tree, levels = trees)) 
  
  baseline = data %>% filter(setting == 'baseline') %>% pull(KS) %>% median
  
  g = ggplot(data, aes(x = level, y = KS, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "KS distance") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
  
  return(g)
}


plot_summary <- function(path, settings, title) {
  
  simStats = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'simStats.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  hpdTrees = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'hpdTreesOutput.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  mccP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'mccOutput.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  g1 = plot_diversity(simStats, settings)
  g2 = plot_wRF_boxplot(mccP, settings)
  g3 = plot_shPI_boxplot(mccP, settings)
  g4 = plot_KS_boxplot(mccP, settings)
  g5 = plot_treecov_lollipop(hpdTrees, settings)

  g = (g1 + ggtitle(title) + theme(plot.title.position = "plot")) + guide_area() +
    g2 + g4 + g3 + g5 + plot_layout(ncol = 2, nrow = 3, guides = 'collect')
  
  return(g)
}


p_td = plot_summary(path_td, settings_td, 'A: non-sequential recordings')
p_tw = plot_summary(path_tw, settings_tw, 'B: sequential recordings')

pdf('pdf/SuppFigure7_ExpVariationsTrees.pdf', height = 10, width = 11)
p_td
p_tw
dev.off()
