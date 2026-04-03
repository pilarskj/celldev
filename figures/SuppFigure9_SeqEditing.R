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

## tree inference
mcc_data = lapply(mccP, function(x) {x %>% 
    bind_rows(.id = 'tree')}) %>%
  bind_rows(.id = 'setting') %>%
  mutate(tree = factor(tree, levels = trees))


p_wRF = ggplot(mcc_data, aes(x = setting, y = wRF, color = tree)) + 
  geom_boxplot(outlier.size = 0.5) + # with legend
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = NULL, y = "Weighted RF distance", color = NULL)

p_RF = ggplot(mcc_data, aes(x = setting, y = RF, color = tree)) + 
  geom_boxplot(outlier.size = 0.5) + # with legend
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = NULL, y = "RF distance", color = NULL)

p_shPI = ggplot(mcc_data, aes(x = setting, y = shPI, color = tree)) + 
  geom_boxplot(outlier.size = 0.5) + # with legend
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = NULL, y = "Shared PI", color = NULL)

p_WS = ggplot(mcc_data, aes(x = setting, y = WS, color = tree)) + 
  geom_boxplot(outlier.size = 0.5) + 
  scale_y_continuous(limits = c(0, 2.5)) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = NULL, y = "Wasserstein distance", color = NULL)

p_KS = ggplot(mcc_data, aes(x = setting, y = KS, color = tree)) + 
  geom_boxplot(show.legend = F, outlier.size = 0.5) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = NULL, y = "KS distance", color = NULL)


# optional: simulation statistics
sim_data = lapply(simStats, function(x) {x %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, barcodeDiv, saturation)}) %>%
  bind_rows(.id = 'setting') %>%
  mutate(tree = factor(tree, levels = trees))
  
p_sim = ggplot(sim_data, aes(x = setting, y = barcodeDiv, color = tree)) + 
  geom_boxplot(show.legend = F, outlier.size = 0.5) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette, labels = tree_labels) +
  labs(x = NULL, y = "Prop. unique barcodes (%)")


## phylodynamic inference
bias_data = lapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
p_bias = ggplot(bias_data, aes(x = setting, y = bias, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.1) + # line at 0
    scale_x_discrete(labels = params_labels) +
    expand_limits(y = c(-0.5, 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias") +
    facet_grid(rows = vars(parameter), labeller = labeller(parameter = params_labels), scales = "free_y")
  

hpd_data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
      mutate(parameter = factor(parameter, levels = params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
p_hpd = ggplot(hpd_data, aes(x = setting, y = hpd_width, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_x_discrete(labels = params_labels) +
    expand_limits(y = c(0, 1)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative HPD width", color = NULL) +
    facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
  

# optional?
prop_params = c('birthRate', 'deathRate')
prop_data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_proportion') %>% 
      bind_rows(.id = 'tree') %>%
      select(tree, seed, all_of(prop_params)) %>%
      pivot_longer(cols = all_of(prop_params), names_to = 'parameter', values_to = 'hpd_prop') %>%
      mutate(parameter = factor(parameter, levels = prop_params),
             tree = factor(tree, levels = trees)) }) %>%
    bind_rows(.id = 'setting') 
  
p_prop = ggplot(prop_data, aes(x = setting, y = hpd_prop, color = tree)) + 
    geom_boxplot(show.legend = F, outlier.size = 0.5) + 
    scale_x_discrete(labels = params_labels) +
    expand_limits(y = 0) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "HPD proportion", color = NULL) +
    facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
  

(((p_wRF + ggtitle('A')) / p_RF / p_shPI / p_WS / p_KS) |
  (p_bias + ggtitle('B')) | p_hpd) +
  plot_layout(guides = 'collect', axes = "collect") & 
  theme(plot.title.position = "plot",
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        legend.position = "bottom", legend.direction = "vertical", legend.justification = "right")
ggsave('pdf/SuppFigure9_SeqEditing.pdf', height = 12, width = 10)

