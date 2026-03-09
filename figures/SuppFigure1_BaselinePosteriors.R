# Posterior distributions

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

# paths
analysis_td = paste0(dir, "/TiDe/baseline/analysisOutput")
analysis_tw = paste0(dir, "/Typewriter/baseline/analysisOutput")

tree_params_f = paste0(dir, "/Trees/tree_params.csv")

# true parameters
treeParams = read.csv(tree_params_f)
treeParams = treeParams %>%
  select(tree, birthRate_true = birthRate, deathRate_true = deathRate,
         treeHeight_true = treeHeight, treeLength_true = treeLength) %>%
  mutate(growthRate_true = birthRate_true - deathRate_true) %>%
  mutate(seed = str_extract(tree, '[0-9]+') %>% as.numeric,
         tree = factor(str_extract(tree, 'tree_[a-z]+'), levels = trees)) 

editRate_true = 0.05


plot_posteriors <- function(analysis_dir, title) {
  
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
    #geom_hline(aes(yintercept = priors$birthRate["median"]), linetype = "dashed", color = "grey") + 
    geom_pointrange(aes(ymin = birthRate_lower, ymax = birthRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_step(aes(y = birthRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, by = 0.2)) +
    labs(x = "Simulation", y = "Posterior division rate", color = NULL)
  
  # death rate
  pD = ggplot(data = infOutputDat, aes(x = ID, y = deathRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$deathRate["lower"], ymax = priors$deathRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    #geom_hline(aes(yintercept = priors$deathRate["median"]), linetype = "dashed", color = "grey") + 
    geom_pointrange(aes(ymin = deathRate_lower, ymax = deathRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_step(aes(y = deathRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.05)) +
    labs(x = "Simulation", y = "Posterior death rate", color = NULL) 
  
  # growth rate
  pG = ggplot(data = infOutputDat, aes(x = ID, y = growthRate_median, color = tree)) +
    geom_pointrange(aes(ymin = growthRate_lower, ymax = growthRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_step(aes(y = growthRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
    labs(x = "Simulation", y = "Posterior growth rate", color = NULL) 
  
  # editing rate
  pE = ggplot(data = infOutputDat, aes(x = ID, y = editRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$editRate["lower"], ymax = priors$editRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    #geom_hline(aes(yintercept = priors$editRate["median"]), linetype = "dashed", color = "grey") + 
    geom_pointrange(aes(ymin = editRate_lower, ymax = editRate_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_segment(aes(x = 0, xend = nrow(infOutputDat), y = editRate_true, yend = editRate_true), color = 'black', alpha = 0.5, linewidth = 0.05) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_y_continuous(limits = c(0, 0.26), breaks = seq(0, 0.25, by = 0.05)) +
    labs(x = "Simulation", y = "Posterior editing rate", color = NULL) 
  
  # tree height
  pTH = ggplot(data = infOutputDat %>% arrange(tree, treeHeight_true) %>% mutate(ID = c(1:nrow(infOutputDat))), 
               aes(x = ID, y = treeHeight_median, color = tree)) +
    geom_pointrange(aes(ymin = treeHeight_lower, ymax = treeHeight_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_line(aes(y = treeHeight_true, group = tree), color = 'black', alpha = 0.5) +
    scale_color_manual(values = palette, labels = tree_labels) +
    ylim(c(0, 41)) + 
    labs(x = "Simulation", y = "Posterior tree height", color = NULL)
  
  # tree length
  pTL = ggplot(data = infOutputDat %>% arrange(tree, treeLength_true) %>% mutate(ID = c(1:nrow(infOutputDat))), 
               aes(x = ID, y = treeLength_median, color = tree)) +
    geom_pointrange(aes(ymin = treeLength_lower, ymax = treeLength_upper), size = 0.1, alpha = 0.6, linewidth = 1) + 
    geom_line(aes(y = treeLength_true, group = tree), color = 'black', alpha = 0.5) + 
    scale_color_manual(values = palette, labels = tree_labels) +
    #ylim(c(0, 3200)) + 
    expand_limits(y = 0) +
    labs(x = "Simulation", y = "Posterior tree length", color = NULL) 
  
  plot = pB + ggtitle(title) + theme(plot.title.position = "plot") + 
    pE + pD + pTH + pG + pTL + 
    plot_layout(guides = 'collect', ncol = 2, nrow = 3) &
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")
  
  return(plot)
}


g_td = plot_posteriors(analysis_td, 'A: non-sequential recordings') 
g_tw = plot_posteriors(analysis_tw, 'B: sequential recordings')

pdf('pdf/SuppFigure1_BaselinePosteriors.pdf', height = 10, width = 9) 
g_td
g_tw
dev.off()
