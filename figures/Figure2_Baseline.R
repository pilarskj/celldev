# Figure 2

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

scenarios = c("non-sequential", "sequential")
paths = c("non-sequential" = paste0(dir, "/TiDe/baseline"), "sequential" = paste0(dir, "/Typewriter/baseline"))


# elements in plot
plot_coverage <- function(infP) {
  data = lapply(infP, '[[', 'coverage') %>% 
    bind_rows(.id = 'tree') %>% 
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'coverage') %>% 
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = params)) %>%
    group_by(tree, parameter) %>% 
    summarise(prop = sum(coverage) * 100 / nSim)
  
  g = ggplot(data, aes(x = parameter, y = prop, color = tree)) +
    geom_point(size = 2, position = position_dodge(width = 0.6)) + 
    geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.1) +
    geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 80%
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Coverage (%)", color = NULL) 

  return(g)
}


plot_bias <- function(infP) {
  data = lapply(infP, '[[', 'bias_rel') %>% 
    bind_rows(.id = 'tree') %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = params))
  
  g = ggplot(data, aes(x = parameter, y = bias, color = tree)) +
    geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = 0, linetype = 'dashed') + # line at 0
    ylim(-1, 2) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias", color = NULL)

  return(g)
}


plot_hpd <- function(infP) {
  par =  c('birthRate', 'deathRate', 'editRate')
  par_labels = params_labels[par]
  data = lapply(infP, '[[', 'hpd_proportion') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, all_of(par)) %>%
    pivot_longer(cols = all_of(par), names_to = 'parameter', values_to = 'hpd_proportion') %>%
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = par))

  g = ggplot(data, aes(x = parameter, y = hpd_proportion, color = tree)) +
    geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = 0.1, linetype = 'dashed') + # threshold at 0.2
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    scale_x_discrete(labels = par_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "HPD proportion", color = NULL)
  
  return(g)
}


plot_wRF <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, wRF)
  
  g = ggplot(data, aes(x = tree, y = wRF, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = 0.2, linetype = 'dashed') + # threshold at 0.2
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    theme(axis.ticks.x = element_blank()) + # remove x axis ticks
    labs(x = NULL, y = "Weighted RF distance", color = NULL)

  return(g)
}


plot_summary <- function(path, title) {
  infP = readRDS(paste0(path, "/analysisOutput/infPerformance.Rdat"))
  mccP = readRDS(paste0(path, "/analysisOutput/mccOutput.Rdat"))
  
  g1 = plot_coverage(infP)
  g2 = plot_bias(infP)
  g3 = plot_hpd(infP)
  g4 = plot_wRF(mccP)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 + g4 +
    plot_layout(ncol = 4, nrow = 1, widths = c(2, 2, 1.3, 1))
  return(g)
}


# plot figures
g_ns = plot_summary(paths['non-sequential'], 'A: non-sequential recordings') 
g_s = plot_summary(paths['sequential'], 'B: sequential recordings')

g_ns / g_s +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))
ggsave('pdf/Figure2_Baseline.pdf', height = 8, width = 11)


# # check correlation with tree size (no)
# mcc_df = mccP %>% bind_rows(.id = "tree")
# tree_params = read.csv(paste0(dir, "/Trees/tree_params.csv")) %>%
#   mutate(seed = as.numeric(str_extract(tree, "[0-9]+")),
#          tree = str_extract(tree, "tree_[a-z]+"))
# mcc_df = mcc_df %>% left_join(tree_params %>% select(tree, seed, nCells), by = c("tree", "seed"))
# ggplot(mcc_df, aes(x = nCells, y = wRF, color = tree)) +
#   geom_point()
# cor.test(mcc_df$nCells, mcc_df$wRF)
