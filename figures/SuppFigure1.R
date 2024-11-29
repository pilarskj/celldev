library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(showtext)
library(HDInterval)
library(patchwork)

# paths
code_dir = "~/Projects/celldev"
analysis_td = "/Volumes/stadler/cEvoUnpublished/2023-Julia-CellDev/Part1/TiDe/baseline/analysisOutput"
analysis_tw = "/Volumes/stadler/cEvoUnpublished/2023-Julia-CellDev/Part1/Typewriter/baseline/analysisOutput"

# settings
palette = c("#E69F00", "#56B4E9", "#009E73","#D55E00", "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman') + 
            theme(axis.text.x = element_blank(), axis.title.x = element_blank()))#, legend.margin = margin(0, 0, 0, 50)))

# trees
trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")
tree_labels = c("tree_s" = "synchronous trees",
                "tree_ss" = "synchronous trees with sampling",
                "tree_sd" = "synchronous trees with cell death",
                "tree_sds" = "synchronous trees with cell death and sampling",
                "tree_bd" = "birth-death-sampling trees")

# true parameters
treeParams = read.csv(file.path(code_dir, 'homog/simulation_trees/tree_params.csv'))
treeParams = treeParams %>%
  select(tree, birthRate_true = birthRate, deathRate_true = deathRate, 
         treeHeight_true = treeHeight, treeLength_true = treeLength) %>%
  mutate(seed = str_extract(tree, '[0-9]+') %>% as.numeric) %>%
  mutate(tree = str_remove(tree, '_[0-9]+'))

editRate_true = 0.05

# priors 
priors = list()
priors$editRate = c(hdi(qlnorm, 0.95, meanlog = -3, sdlog = 1), 
                    "median" = qlnorm(0.5, meanlog = -3, sdlog = 1))
priors$birthRate = c(hdi(qlnorm, 0.95, meanlog = -1.5, sdlog = 1),
                     "median" = qlnorm(0.5, meanlog = -1.5, sdlog = 1))
priors$deathRate = c(hdi(qexp, 0.95, rate = 25),
                     "median" = qexp(0.5, rate = 25))

plot_estimates <- function(analysis_dir, title) {
  
  # load posterior parameter estimates
  infOutputDat = readRDS(file.path(analysis_dir, "infOutput.Rdat")) %>%
    bind_rows(.id = "tree") %>% 
    left_join(treeParams, by = c('tree', 'seed')) %>% # join with true params
    mutate(tree = factor(tree, levels = trees))
  # assign ID
  infOutputDat$ID = seq(1:nrow(infOutputDat))
  
  # birth rate
  pB = ggplot(data = infOutputDat, aes(x = ID, y = birthRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$birthRate["lower"], ymax = priors$birthRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    geom_pointrange(aes(ymin = birthRate_lower, ymax = birthRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_step(aes(y = birthRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, by = 0.2)) +
    labs(x = "Simulation", y = "Posterior birth rate", color = "")
  
  # death rate
  pD = ggplot(data = infOutputDat, aes(x = ID, y = deathRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$deathRate["lower"], ymax = priors$deathRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    geom_pointrange(aes(ymin = deathRate_lower, ymax = deathRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_step(aes(y = deathRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.04)) +
    labs(x = "Simulation", y = "Posterior death rate", color = "") 
  
  # editing rate
  pE = ggplot(data = infOutputDat, aes(x = ID, y = editRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$editRate["lower"], ymax = priors$editRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    #geom_pointrange(x = -2, ymin = priors$editRate['lower'], y = priors$editRate['median'], ymax = priors$editRate['upper'], 
    #                color = 'darkgrey', size = 0.4, alpha = 0.5, linewidth = 1) + 
    geom_pointrange(aes(ymin = editRate_lower, ymax = editRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_segment(aes(x = 0, xend = 100, y = editRate_true, yend = editRate_true), color = 'black', alpha = 0.5, linewidth = 0.05) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 0.26), breaks = seq(0, 0.25, by = 0.05)) +
    labs(x = "Simulation", y = "Posterior editing rate", color = "") 
  
  # tree height
  pTH = ggplot(data = infOutputDat %>% arrange(tree, treeHeight_true) %>% mutate(ID = c(1:nrow(infOutputDat))), 
               aes(x = ID, y = treeHeight_median, color = tree)) +
    geom_pointrange(aes(ymin = treeHeight_lower, ymax = treeHeight_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_line(aes(y = treeHeight_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    ylim(c(0, 41)) + 
    labs(x = "Simulation", y = "Posterior tree height", color = "")
  
  # tree length
  pTL = ggplot(data = infOutputDat %>% arrange(tree, treeLength_true) %>% mutate(ID = c(1:nrow(infOutputDat))), 
               aes(x = ID, y = treeLength_median, color = tree)) +
    geom_pointrange(aes(ymin = treeLength_lower, ymax = treeLength_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_line(aes(y = treeLength_true, group = tree), color = 'black', alpha = 0.5) + 
    #geom_point(aes(y = treeLength_true), color = 'black', alpha = 0.5, shape = '-', size = 5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    ylim(c(0, 3200)) + 
    labs(x = "Simulation", y = "Posterior tree length", color = "") 
  
  # # alternative:
  # pTL = ggplot(data = infOutputDat, aes(x = treeLength_true, y = treeLength_median, color = tree)) +
  #   geom_pointrange(aes(ymin = treeLength_lower, ymax = treeLength_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
  #   geom_abline(slope = 1, intercept = 0, alpha = 0.5, color = 'black') +
  #   scale_color_manual(values = palette, labels = tree_labels) +
  #   xlim(c(0, 3200)) + 
  #   ylim(c(0, 3200)) + 
  #   labs(x = "True tree length", y = "Inferred tree length", color = "") 
  
  plot = (pB + ggtitle(title) + theme(plot.title = element_text(face = "bold"), plot.title.position = "plot") + pD) / 
    (pE + pTH + pTL) + 
    plot_layout(guides = 'collect', widths = c(1.5,1.5,1,1,1))
  #plot = pB + pD + guide_area() + pE + pTH + pTL + 
  #  plot_layout(design = "abc\ndef", guides = 'collect', axis_titles = 'collect_x')
  
  return(plot)
}

g_td = plot_estimates(analysis_td, 'A: non-sequential recordings') 
g_tw = plot_estimates(analysis_tw, 'B: sequential recordings')

pdf('figures/SuppFigure1.pdf', height = 18, width = 12)#, units = "in", res = 300)
wrap_plots(g_td, g_tw, ncol = 1, heights = c(1, 1)) & 
  theme(legend.position = 'bottom', 
        legend.direction = "vertical", 
        legend.justification = "right")
dev.off()

