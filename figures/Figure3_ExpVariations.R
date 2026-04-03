# Figure 3

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

path_td = paste0(dir, "/TiDe")
path_tw = paste0(dir, "/Typewriter")


plot_wRF <- function(path, settings) {
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
    geom_boxplot(outlier.size = 0.5) + 
    geom_hline(yintercept = baseline, linetype = 'dashed') +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Weighted RF distance") +
    facet_grid(cols = vars(manipulation), scales = 'free_x', space = 'free_x')
    
  return(g)
}


plot_hpd <- function(path, settings, variation, baseline_level) {
  
  palette_lev = c("#D55E00","#E69F00", "#009E73")
  
  infP = sapply(settings$setting, function(s) 
    readRDS(file.path(path, s, 'analysisOutput', 'infPerformance.Rdat')), 
    simplify = F, USE.NAMES = T)
  
  data = sapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params))}, simplify = F) %>%
    bind_rows(.id = 'setting') %>%
    left_join(settings) 
  
  # just plot one variation
  data = data %>%
    filter(manipulation %in% c('baseline', variation)) %>%
    mutate(level = ifelse(manipulation == 'baseline', baseline_level, as.numeric(level))) %>%
    mutate(level = as.factor(level))
  
  g = ggplot(data, aes(x = parameter, y = hpd_width, color = level)) +
    geom_boxplot(outlier.size = 0.5) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette_lev) +
    labs(x = NULL, y = "Relative HPD width", color = variation) +
    theme(legend.position = "inside", legend.position.inside = c(0.8, 0.9),
          legend.background = element_rect(linetype = "dashed", color = 'black', linewidth = 0.2),
          legend.title = element_text(size = 10), legend.text = element_text(size = 8),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), axis.text.y = element_text(size = 8))
  
  return(g)
}


# comparison across settings per recorder
pA1 = plot_wRF(path_td, settings_td)
pA2 = plot_hpd(path_td, settings_td, 'editing rate', 0.05)
pA3 = plot_hpd(path_td, settings_td, 'no. targets', 20)
pB1 = plot_wRF(path_tw, settings_tw)
pB2 = plot_hpd(path_tw, settings_tw, 'editing rate', 0.05)
pB3 = plot_hpd(path_tw, settings_tw, 'no. tapes', 20)


p_td = pA1 + ggtitle('A: non-sequential recordings') + theme(plot.title.position = "plot") + pA2 + pA3 + plot_layout(ncol = 3, widths = c(1.75, 1, 1), axis_titles = "collect")
p_tw = pB1 + ggtitle('B: sequential recordings') + theme(plot.title.position = "plot") + pB2 + pB3 + plot_layout(ncol = 3, widths = c(1.75, 1, 1), axis_titles = "collect")
p = p_td / p_tw


# associations across ALL settings
paths_td = sapply(settings_td$setting, function(s) file.path(path_td, s))
names(paths_td) = paste0("TiDe_", names(paths_td))
paths_tw = sapply(settings_tw$setting, function(s) file.path(path_tw, s))
names(paths_tw) = paste0("Typewriter_", names(paths_tw))
paths = c(paths_td, paths_tw)

# recording statistics
simStats = sapply(paths, function(p) 
  readRDS(file.path(p, 'analysisOutput', 'simStats.Rdat')), 
  simplify = F, USE.NAMES = T)

dfsim = sapply(simStats, function(x) {x %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, barcodeDiv, nEdited)}, simplify = F) %>%
  bind_rows(.id = 'setting') %>%
  mutate(propUnique = barcodeDiv * 100)

# parameter inference
infP = sapply(paths, function(p) 
  readRDS(file.path(p, 'analysisOutput', 'infPerformance.Rdat')), 
  simplify = F, USE.NAMES = T)

# extract metrics for growth rate
dfbias = sapply(infP, function(x) {
  lapply(x, '[[', 'bias_rel') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, bias = growthRate)}, simplify = F) %>%
  bind_rows(.id = 'setting') 
dfhpd = sapply(infP, function(x) {
  lapply(x, '[[', 'hpd_width_rel') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, hpd = growthRate)}, simplify = F) %>%
  bind_rows(.id = 'setting') 

# tree inference
mccP = sapply(paths, function(p) 
  readRDS(file.path(p, 'analysisOutput', 'mccOutput.Rdat')), 
  simplify = F, USE.NAMES = T)

dfmcc = sapply(mccP, function(x) {x %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, wRF, shPI, KS)}, simplify = F) %>%
  bind_rows(.id = 'setting') %>% 
  mutate(seed = as.character(seed))

df = dfsim %>% full_join(dfbias) %>% full_join(dfhpd) %>% full_join(dfmcc) %>%
  mutate(tree = factor(tree, levels = trees)) 

g1 = ggplot(df, aes(x = propUnique, y = wRF)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,100), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Prop. unique barcodes (%)", y = "Weighted RF distance", color = NULL)

g2 = ggplot(df, aes(x = wRF, y = bias)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(-0.5,0.5), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Weighted RF distance", y = "Relative bias\n(growth rate)", color = NULL)

g3 = ggplot(df, aes(x = wRF, y = hpd)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,1), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0,1.5), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Weighted RF distance", y = "Relative HPD width\n(growth rate)", color = NULL)

g = g1 + ggtitle("C") + theme(plot.title.position = "plot") + g2 + g3 + plot_layout(guides = "collect") &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")

# final plot
p / g
ggsave('pdf/Figure3_ExpVariations.pdf', height = 12, width = 11) 


### add to supplement:
gPI = ggplot(df, aes(x = propUnique, y = shPI)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,100), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0,1.1), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Prop. unique barcodes (%)", y = "Shared PI", color = NULL)

gKS = ggplot(df, aes(x = propUnique, y = KS)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,100), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(0,1.1), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Prop. unique barcodes (%)", y = "KS distance", color = NULL)

# biases to division and death 
dfbias = sapply(infP, function(x) {
  lapply(x, '[[', 'bias_rel') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, birthRate, deathRate)}, simplify = F) %>%
  bind_rows(.id = 'setting') 
df = dfsim %>% full_join(dfbias) %>% full_join(dfmcc) %>%
  mutate(tree = factor(tree, levels = trees)) 

gBiasDiv = ggplot(df, aes(x = propUnique, y = birthRate)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,100), breaks = pretty_breaks()) +
  scale_y_continuous(limits = c(-0.5,0.5), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Prop. unique barcodes (%)", y = "Relative bias\n(division rate)", color = NULL)

gBiasDeath = ggplot(df, aes(x = propUnique, y = deathRate)) +
  geom_point(aes(color = tree), size = 0.5) +
  stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.2, color = "black") +
  stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  scale_x_continuous(limits = c(0,100), breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = "Prop. unique barcodes (%)", y = "Relative bias\n(death rate)", color = NULL)

(gPI | gKS) / (gBiasDiv | gBiasDeath) + plot_layout(guides = "collect") &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")

ggsave('pdf/SuppFigure8_DiversityBiases.pdf', height = 8, width = 8) 
