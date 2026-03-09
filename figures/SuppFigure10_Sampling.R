# Inference performance on large vs. subsampled trees

setwd("~/Projects/celldev/figures")
source("homog_functions.R")
theme_update(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
             plot.title.position = "plot",
             legend.position = "bottom", legend.direction = "vertical", legend.justification = "right")

# paths
paths_td = c("large" = paste0(dir, "/TiDe/baseline_large/analysisOutput"), 
             "sampled" = paste0(dir, "/TiDe/baseline_sampled/analysisOutput"))

paths_tw = c("large" = paste0(dir, "/Typewriter/baseline_large/analysisOutput"), 
             "sampled" = paste0(dir, "/Typewriter/baseline_sampled/analysisOutput"))

settings = c("large", "sampled")
setting_labels = c("large" = "larger trees\n(10% sampling)", 
                   "sampled" = "smaller trees\n(1% sampling)")

# only look at trees with sampling
nRuns = 60
palette = c("#56B4E9", "#D55E00", "#CC79A7")
trees = c("tree_ss", "tree_sds", "tree_bd")
tree_labels = c("tree_ss" = "synchronous trees",
                "tree_sds" = "synchronous trees with cell death",
                "tree_bd" = "birth-death trees")


# tree inference
plot_wRF <- function(data) {
  
  g = ggplot(data, aes(x = setting, y = wRF, color = tree)) + 
    geom_boxplot(outlier.size = 0.5) + 
    scale_x_discrete(labels = setting_labels) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Weighted RF distance", color = NULL)
  
  return(g)
}


plot_shPI <- function(data) {
  
  g = ggplot(data, aes(x = setting, y = shPI, color = tree)) + 
    geom_boxplot(outlier.size = 0.5) + 
    scale_x_discrete(labels = setting_labels) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared PI", color = NULL)
  
  return(g)
}


plot_KS <- function(data) {
  
  g = ggplot(data, aes(x = setting, y = KS, color = tree)) + 
    geom_boxplot(outlier.size = 0.5) + 
    scale_x_discrete(labels = setting_labels) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "KS distance", color = NULL)
  
  return(g)
}


# parameter inference
plot_bias <- function(infP, converged_runs) {
  
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'bias_rel') %>% 
      bind_rows(.id = 'tree') %>%
      filter(tree %in% trees) %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') 
  }) %>% bind_rows(.id = 'setting') %>%
    mutate(parameter = factor(parameter, levels = params),
           tree = factor(tree, levels = trees),
           setting = factor(setting, levels = settings),
           seed = as.numeric(seed)) %>%
    left_join(converged_runs) %>% 
    filter(converged)
  
  g = ggplot(data, aes(x = parameter, y = bias, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    expand_limits(y = c(-0.2, 0.2)) +
    scale_x_discrete(labels = params_labels) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias", color = NULL) +
    facet_wrap(~setting, nrow = 1, labeller = as_labeller(setting_labels)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))
    #facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
 
  return(g) 
}


plot_hpd <- function(infP, converged_runs) {
  
  data = lapply(infP, function(x) {
    lapply(x, '[[', 'hpd_width_rel') %>% 
      bind_rows(.id = 'tree') %>%
      filter(tree %in% trees) %>%
      select(tree, seed, all_of(params)) %>%
      pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'hpd_width') 
  }) %>% bind_rows(.id = 'setting') %>%
    mutate(parameter = factor(parameter, levels = params),
           tree = factor(tree, levels = trees),
           setting = factor(setting, levels = settings),
           seed = as.numeric(seed)) %>%
    left_join(converged_runs) %>% 
    filter(converged)
  
  g = ggplot(data, aes(x = parameter, y = hpd_width, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    expand_limits(y = c(0, 0.5)) +
    scale_x_discrete(labels = params_labels) +
    scale_y_continuous(breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative HPD width", color = NULL) +
    facet_wrap(~setting, nrow = 1, labeller = as_labeller(setting_labels)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))
    #facet_grid(rows = vars(parameter), scales = 'free_y', labeller = labeller(parameter = params_labels))
  
  return(g) 
}


plot_summary <- function(paths, title) {
  
  # remove runs which did not converge
  converged_runs = lapply(paths, function(p) {
    dat = readRDS(file.path(p, "infOutput.Rdat")) %>%
      bind_rows(.id = "tree") %>% 
      filter(tree %in% trees) %>%
      mutate(tree = factor(tree, levels = trees)) %>%
      mutate(converged = ifelse(minESS >= 200, T, F)) %>%
      select(tree, seed, converged)
    return(dat)
  }) %>% bind_rows(.id = 'setting')
  
  mccP = sapply(paths, function(p) readRDS(file.path(p, 'mccOutput.Rdat')), simplify = F, USE.NAMES = T)
  infP = sapply(paths, function(p) readRDS(file.path(p, 'infPerformance.Rdat')), simplify = F, USE.NAMES = T)
  
  mcc_dat = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') }) %>%
    bind_rows(.id = 'setting') %>%
    filter(tree %in% trees) %>%
    mutate(tree = factor(tree, levels = trees),
           setting = factor(setting, levels = settings)) %>%
    left_join(converged_runs) %>% 
    filter(converged)
  
  p1 = plot_wRF(mcc_dat)
  p2 = plot_shPI(mcc_dat)
  p3 = plot_KS(mcc_dat)
  p4 = plot_bias(infP, converged_runs)
  p5 = plot_hpd(infP, converged_runs)
  
  #leg <<- get_legend(p1)
  
  plot = (p1 + ggtitle(title) + p2 + p3) /
    (p4 + p5) +
    plot_layout(guides = 'collect') #& theme(legend.position = "none")
  
  return(plot)
}


p_td = plot_summary(paths_td, 'A: non-sequential recordings')
p_tw = plot_summary(paths_tw, 'B: sequential recordings')

pdf('pdf/SuppFigure10_Sampling.pdf', height = 8, width = 9) 
p_td
p_tw
dev.off()


