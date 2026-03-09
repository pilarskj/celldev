# Figure 4

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

# update plotting options
palette = c("non-sequential" = "#56B4E9", "sequential" = "#CC79A7")

# data
paths = c('non-sequential' = file.path(dir, 'TiDe/baseline/analysisOutput'),
          'sequential' = file.path(dir, 'Typewriter/sequential_betterPrior/analysisOutput'))

mccP = sapply(paths, function(p) readRDS(file.path(p, 'mccOutput.Rdat')), simplify = F, USE.NAMES = T)
infP = sapply(paths, function(p) readRDS(file.path(p, 'infPerformance.Rdat')), simplify = F, USE.NAMES = T)


## tree inference (MCC)
mcc_data = lapply(mccP, function(x) {x %>% 
    bind_rows(.id = 'tree') }) %>%
  bind_rows(.id = 'setting') 

# topology and branches
g1 = ggplot(mcc_data, aes(x = setting, y = wRF)) + 
  geom_boxplot(aes(color = setting), outlier.size = 0.5) + 
  geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 'wilcox.test', test.args = list(paired = T, alternative = "greater"), map_signif_level = T) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "Weighted RF distance", tag = "A")
  
# topology only
g2 = ggplot(mcc_data, aes(x = setting, y = shPI)) + 
  geom_boxplot(aes(color = setting), outlier.size = 0.5) + 
  geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 'wilcox.test', test.args = list(paired = T, alternative = "less"), map_signif_level = T) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "Shared PI")

# branches only
g3 = ggplot(mcc_data, aes(x = setting, y = KS)) + 
  geom_boxplot(aes(color = setting), outlier.size = 0.5) + 
  geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 'wilcox.test', test.args = list(paired = T, alternative = "greater"), map_signif_level = T) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "KS distance")

# g1 + g2 + g3 + plot_layout(guides = "collect") & theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) 


## phylodynamic inference
params = c('birthRate', 'deathRate', 'growthRate')
params_labels = params_labels[params]

bias = lapply(infP, function(x) {
  lapply(x, '[[', 'bias_rel') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, all_of(params)) %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias_rel') %>%
    mutate(parameter = factor(parameter, levels = params)) }) %>%
  bind_rows(.id = 'setting') 

gBias = ggplot(bias, aes(x = setting, y = bias_rel)) + 
  geom_boxplot(aes(color = setting), outlier.size = 0.5) + 
  geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 'wilcox.test', test.args = list(paired = T, alternative = "greater"), map_signif_level = T) +
  scale_y_continuous(limits = c(-1, 2.2), breaks = scales::pretty_breaks()) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "Relative bias", tag = "B") +
  facet_wrap(vars(parameter), labeller = as_labeller(params_labels))

hpd = lapply(infP, function(x) {
  lapply(x, '[[', 'hpd_width_rel') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, all_of(params)) %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') %>%
    mutate(parameter = factor(parameter, levels = params)) }) %>%
  bind_rows(.id = 'setting') 

gHPD = ggplot(hpd, aes(x = setting, y = hpd_width)) + 
  geom_boxplot(aes(color = setting), outlier.size = 0.5) + 
  geom_signif(comparisons = list(c("non-sequential", "sequential")), test = 'wilcox.test', test.args = list(paired = T, alternative = "greater"), map_signif_level = T) +
  scale_y_continuous(limits = c(0, 8.5), breaks = scales::pretty_breaks()) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "Relative HPD width") +
  facet_wrap(vars(parameter), labeller = as_labeller(params_labels))

# Wilcoxon test
# trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")
# for (t in trees) {
#    print(t)
#    ns = hpd %>% filter(tree == t, parameter == "birthRate", setting == "non-sequential") %>% pull(hpd_width)
#    s = hpd %>% filter(tree == t, parameter == "birthRate", setting == "sequential") %>% pull(hpd_width)
#    print(wilcox.test(ns, s, paired = TRUE, alternative = "greater"))
# }
ns = mcc_data %>% filter(setting == "non-sequential") %>% pull(wRF)
s = mcc_data %>% filter(setting == "sequential") %>% pull(wRF)
#ns = hpd %>% filter(parameter == "deathRate", setting == "non-sequential") %>% pull(hpd_proportion)
#s = hpd %>% filter(parameter == "deathRate", setting == "sequential") %>% pull(hpd_proportion)
wilcox.test(ns, s, paired = TRUE, alternative = "greater")

(g1 + g2 + g3) / (gBias + gHPD) + plot_layout(guides = "collect") & 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom", legend.title = element_blank()) 
ggsave('pdf/Figure4_SeqEditing.pdf', height = 6, width = 8)


