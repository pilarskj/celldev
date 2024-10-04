library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(showtext)
library(HDInterval)
font_add("lmroman", regular = "~/Projects/celldev/Part1/Results/lmroman10-regular.otf")

trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")

dir = '~/Projects/celldev_data/analysisOutput'
inf_f = file.path(dir, 'infOutput.Rdat')

# create table and filter valid simulations
infOutput = readRDS(inf_f)
infOutputDat = bind_rows(infOutput, .id = "tree")
infOutputDat = infOutputDat %>% 
  mutate(tree = factor(tree, levels = trees)) %>% 
  filter(minESS > 200)
infOutputDat %>% group_by(tree) %>% summarise(validSim = n()) 

bRout = infOutputDat %>% 
  select(tree, seed, starts_with('birthRate')) %>%
  group_by(tree) %>%
  summarise(across(starts_with('birthRate'), ~mean(.x)))

ggplot(bRout, aes(x = tree, y = birthRate_median)) + 
  geom_pointrange(aes(ymin = birthRate_lower, ymax = birthRate_upper), size = 0.4) 


perf_f = file.path(dir, 'infPerformance.Rdat')
infP = readRDS(perf_f)

bias_rel = lapply(infP, '[[', 'bias_rel') %>%  bind_rows(.id = "tree")

ggplot(bias_rel, aes(x = tree, y = birthRate)) + 
  geom_boxplot() +
  ylim(c(-1,1)) +
  geom_hline(yintercept = 0, linetype = "dotted") 


bias_rel = bias_rel %>% 
  select(tree, seed, birthRate) %>%
  group_by(tree) %>% 
  summarise(mean = mean(birthRate), sd = sd(birthRate))

p = ggplot(bias_rel, aes(x = tree, y = mean)) + 
  geom_pointrange(aes(ymin = mean - sd/2, ymax = mean + sd/2), size = 0.4) +
  ylim(c(-1,1)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw(base_family = "lmroman", base_size = 16) #+ 


infOutputDat$ID = seq(1:nrow(infOutputDat))
  
# Test squeezed figure
ggplot(data = infOutputDat, aes(x = ID, y = editRate_median, color = tree)) +
  geom_pointrange(aes(ymin = editRate_lower, ymax = editRate_upper), size = 0.1, alpha = 0.4, linewidth = 1) + 
  geom_hline(yintercept = 0.45) +  
  ylim(c(0, 1)) + 
  # theme_bw(base_family = "lmroman", base_size = 20) + 
  xlab("Simulation") + 
  ylab("Posterior editing rate")

phyloParams = left_join(infOutputDat, phyloParams, by = 'tree')

ggplot(data = infOutputDat, aes(x = ID, y = deathRate_median, color = tree)) +
  geom_pointrange(aes(ymin = deathRate_lower, ymax = deathRate_upper), size = 0.1, alpha = 0.5, linewidth = 1) + 
  geom_pointrange(x = -2, ymin = priors$deathRate['lower'], y = priors$deathRate['median'], ymax = priors$deathRate['upper'], 
                  color = 'darkgrey', size = 0.2, alpha = 0.5, linewidth = 1) + 
  geom_step(data = phyloParams, aes(x = ID, y = deathRate), color = 'black', alpha = 0.5) +
  ylim(c(0, 0.2)) + 
  xlim(c(-2, 100)) +
  # theme_bw(base_family = "lmroman", base_size = 20) + 
  xlab("Simulation") + 
  ylab("Posterior death rate")

geom_hline(data = phyloParams, aes(yintercept = birthRate)) +

ggplot(data = infOutputDat, aes(x = ID, y = treeHeight_median, color = tree)) +
  geom_pointrange(aes(ymin = treeHeight_lower, ymax = treeHeight_upper), size = 0.1, alpha = 0.4, linewidth = 1) + 
  geom_hline(yintercept = 0.45) +  
  ylim(c(0, 1)) + 
  # theme_bw(base_family = "lmroman", base_size = 20) + 
  xlab("Simulation") + 
  ylab("Posterior editing rate")


palette <- c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
params = c('editRate', 'treeHeight', 'treeLength', 'birthRate', 'deathRate')
trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")

# lollipop chart for plotting coverage
cov = lapply(infPerformance, '[[', 'coverage') %>% 
  bind_rows(.id = 'tree') %>% 
  pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'coverage') %>% 
    group_by(tree, parameter) %>% 
  summarise(prop = sum(coverage) / nSim)
ggplot(data = cov, aes(x = parameter, y = prop, color = tree)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.2) +
  scale_color_manual(values = palette) +
  geom_hline(yintercept = 0.8, linetype = 'dashed') 


bias = lapply(infPerformance, '[[', 'bias_rel') %>% 
  bind_rows(.id = 'tree') %>%
  mutate(tree = factor(tree, levels = trees)) %>%
  pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias')
#bias = bias %>% group_by(tree, parameter) %>% summarise(mean = mean(bias), sd = sd(bias))
ggplot(data = bias, aes(x = parameter, y = bias, color = tree)) +
  geom_violin(position = position_dodge(0.6)) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 0.5), geom = "pointrange", size = 0.1, position = position_dodge(0.6)) +
  scale_color_manual(values = palette) +
  geom_hline(yintercept = 0, linetype = 'dashed') 

hpd_width = lapply(infPerformance, '[[', 'hpd_width_rel') %>% 
  bind_rows(.id = 'tree') %>%
  mutate(tree = factor(tree, levels = trees)) %>%
  pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width')
#hpd_width = hpd_width %>% group_by(tree, parameter) %>% summarise(mean = mean(hpd_width), sd = sd(hpd_width))
ggplot(data = hpd_width, aes(x = parameter, y = hpd_width, color = tree)) +
  geom_violin(position = position_dodge(0.6)) +
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 0.5), geom = "pointrange", size = 0.1, position = position_dodge(0.6)) +
  scale_color_manual(values = palette) +
  geom_hline(yintercept = 0, linetype = 'dashed') 



## test figures for mcc
mcc_f = file.path(dir, 'mccOutput.Rdat')
mccOutput = readRDS(mcc_f)
metrics = c('wRF', 'Nye', 'shPI')
mcc = mccOutput %>% 
  bind_rows(.id = 'tree') %>%
  mutate(tree = factor(tree, levels = trees)) 
  #pivot_longer(cols = all_of(metrics), names_to = 'metric', values_to = 'value')
  select(tree, seed, wRF)
ggplot(data = mcc, aes(y = wRF, color = tree)) + # x = metric, y = value
  geom_boxplot() + # position = position_dodge(0.6)
  scale_color_manual(values = palette) +
  ylim(0, 1) +
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()) +
  geom_hline(yintercept = 0.2, linetype = 'dashed') 
