library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(showtext)
library(HDInterval)
font_add("lmroman", regular = "~/Desktop/Master Thesis/Code and Data/Part1/Results/lmroman10-regular.otf")

trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")

out_f = '~/Desktop/Master Thesis/Code and Data/Part1/Typewriter/baseline/analysisOutput/infOutput.Rdat'

# create table and filter valid simulations
infOutput = readRDS(out_f)
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


perf_f = '~/Desktop/Master Thesis/Code and Data/Part1/Typewriter/baseline/analysisOutput/infPerformance.Rdat'
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


