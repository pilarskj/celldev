# functions for all figures

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(HDInterval)
library(ggplot2)
library(showtext)
library(patchwork)
library(cowplot)
library(ggsignif)
library(ggpubr)
library(purrr)
library(ggridges)
library(scales)

# plotting options
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman'))

# data directory
dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-CellDev/homog"

# trees
trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")
tree_labels = c("tree_s" = "synchronous trees",
                "tree_ss" = "synchronous trees with sampling",
                "tree_sd" = "synchronous trees with cell death",
                "tree_sds" = "synchronous trees with cell death and sampling",
                "tree_bd" = "birth-death trees with sampling")

# number of simulations per tree
nSim = 20

# priors 
priors = list()
priors$editRate = c(hdi(qlnorm, 0.95, meanlog = -3, sdlog = 1), 
                    "median" = qlnorm(0.5, meanlog = -3, sdlog = 1))
priors$birthRate = c(hdi(qlnorm, 0.95, meanlog = -1.5, sdlog = 1),
                     "median" = qlnorm(0.5, meanlog = -1.5, sdlog = 1))
priors$deathRate = c(hdi(qexp, 0.95, rate = 25),
                     "median" = qexp(0.5, rate = 25))

# parameters 
params = c('birthRate', 'deathRate', 'growthRate', 'editRate', 'treeHeight', 'treeLength') 
params_labels = c('birthRate' = 'division rate', 
                  'deathRate' = 'death rate',
                  'growthRate' = 'growth rate',
                  'editRate' = 'editing rate',
                  'treeHeight' = 'tree height', 
                  'treeLength' = 'tree length')


