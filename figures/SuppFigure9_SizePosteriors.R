# Compare posteriors at different tree sizes

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

# paths
paths_td = c("large" = paste0(dir, "/TiDe/baseline_large/analysisOutput"), 
             "sampled" = paste0(dir, "/TiDe/baseline_sampled/analysisOutput"))

paths_tw = c("large" = paste0(dir, "/Typewriter/baseline_large/analysisOutput"), 
             "sampled" = paste0(dir, "/Typewriter/baseline_sampled/analysisOutput"))

tree_params_fs = c("large" =  paste0(dir, "/LargeTrees/tree_params.csv"),
                   "sampled" = paste0(dir, "/SampledTrees/tree_params.csv"))

setting_labels = c("large" = "larger trees (10% sampling)", 
                   "sampled" = "smaller trees (1% sampling)")

# only look at trees with sampling
nRuns = 60
palette = c("#56B4E9", "#D55E00", "#CC79A7")
trees = c("tree_ss", "tree_sds", "tree_bd")
tree_labels = c("tree_ss" = "synchronous trees",
                "tree_sds" = "synchronous trees with cell death",
                "tree_bd" = "birth-death trees")

# only show phylodynamic parameters and editing rate
params = c("birthRate", "deathRate", "growthRate", "editRate")
params_labels = params_labels[params]

# get true parameters
treeParams = lapply(tree_params_fs, function(f) {
  # true parameters
  df = read.csv(f)
  df = df %>%
    select(tree, birthRate_true = birthRate, deathRate_true = deathRate) %>%
    mutate(growthRate_true = birthRate_true - deathRate_true) %>%
    mutate(seed = str_extract(tree, '[0-9]+') %>% as.numeric,
           tree = str_extract(tree, 'tree_[a-z]+')) %>%
    filter(tree %in% trees) %>%
    mutate(tree = factor(tree, levels = trees))
  return(df)
})  
# should be the same!
all(treeParams[["large"]] == treeParams[["sampled"]])
treeParams = treeParams[["large"]]

editRate_true = 0.05


plot_posteriors <- function(paths, title) {
  
  # load posterior parameter estimates
  infOutputDat = lapply(paths, function(p) {
    dat = readRDS(file.path(p, "infOutput.Rdat")) %>%
      bind_rows(.id = "tree") %>% 
      filter(tree %in% trees) %>%
      left_join(treeParams, by = c('tree', 'seed')) %>% # join with true params
      mutate(tree = factor(tree, levels = trees)) %>%
      mutate(converged = ifelse(minESS >= 200, T, F))
    # assign ID
    dat$ID = seq(1:nrow(dat))
    return(dat)
  }) %>% bind_rows(.id = 'setting')
  
  # birth rate
  pB = ggplot(data = infOutputDat, aes(x = ID, y = birthRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$birthRate["lower"], ymax = priors$birthRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    #geom_hline(aes(yintercept = priors$birthRate["median"]), linetype = "dashed", color = "grey") + 
    geom_pointrange(aes(ymin = birthRate_lower, ymax = birthRate_upper, alpha = converged), size = 0.1, linewidth = 1) +
    geom_step(aes(y = birthRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_y_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, by = 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 0.1)) +
    labs(x = "Simulation", y = "Posterior division rate", color = NULL) +
    facet_wrap(~setting, labeller = as_labeller(setting_labels))
  
  # death rate
  pD = ggplot(data = infOutputDat, aes(x = ID, y = deathRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$deathRate["lower"], ymax = priors$deathRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    #geom_hline(aes(yintercept = priors$deathRate["median"]), linetype = "dashed", color = "grey") + 
    geom_pointrange(aes(ymin = deathRate_lower, ymax = deathRate_upper, alpha = converged), size = 0.1, linewidth = 1) + 
    geom_step(aes(y = deathRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_y_continuous(limits = c(0, 0.2), breaks = seq(0, 0.2, by = 0.05)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 0.1)) +
    labs(x = "Simulation", y = "Posterior death rate", color = NULL) +
    facet_wrap(~setting) +
    theme(strip.background = element_blank(), strip.text = element_blank())
  
  # growth rate
  pG = ggplot(data = infOutputDat, aes(x = ID, y = growthRate_median, color = tree)) +
    geom_pointrange(aes(ymin = growthRate_lower, ymax = growthRate_upper, alpha = converged), size = 0.1, linewidth = 1) + 
    geom_step(aes(y = growthRate_true, group = tree), color = 'black', alpha = 0.5) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, by = 0.1)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 0.1)) +
    labs(x = "Simulation", y = "Posterior growth rate", color = NULL) +
    facet_wrap(~setting) +
    theme(strip.background = element_blank(), strip.text = element_blank())
  
  # editing rate
  pE = ggplot(data = infOutputDat, aes(x = ID, y = editRate_median, color = tree)) +
    geom_ribbon(aes(ymin = priors$editRate["lower"], ymax = priors$editRate["upper"]), 
                alpha = 0.2, fill = "grey", color = 'grey', linewidth = 0) +
    #geom_hline(aes(yintercept = priors$editRate["median"]), linetype = "dashed", color = "grey") + 
    geom_pointrange(aes(ymin = editRate_lower, ymax = editRate_upper, alpha = converged), size = 0.1, linewidth = 1) + 
    geom_segment(aes(x = 0, xend = nRuns, y = editRate_true, yend = editRate_true), color = 'black', alpha = 0.5, linewidth = 0.05) +
    scale_y_continuous(limits = c(0, 0.26), breaks = seq(0, 0.25, by = 0.05)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    scale_alpha_manual(values = c("TRUE" = 0.6, "FALSE" = 0.1)) +
    labs(x = "Simulation", y = "Posterior editing rate", color = NULL) +
    facet_wrap(~setting) +
    theme(strip.background = element_blank(), strip.text = element_blank()) 

  plot = (pB + ggtitle(title) + theme(plot.title.position = "plot")) + pD + pG + pE + 
            plot_layout(guides = 'collect', ncol = 1) &
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 8))

  return(plot)
}

p_td = plot_posteriors(paths_td, 'A: non-sequential recordings')
p_tw = plot_posteriors(paths_tw, 'B: sequential recordings')

pdf('pdf/SuppFigure9_SizePosteriors.pdf', height = 9, width = 9) 
p_td
p_tw
dev.off()
