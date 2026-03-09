# Additional tree metrics at baseline

setwd("~/Projects/celldev/figures")
source("homog_functions.R")
theme_update(axis.ticks.x = element_blank(), axis.text.y = element_text(size = 8))

scenarios = c("non-sequential", "sequential")
paths = c("non-sequential" = paste0(dir, "/TiDe/baseline/analysisOutput"), "sequential" = paste0(dir, "/Typewriter/baseline/analysisOutput"))


# plots for HPD trees
plot_coverage_lollipop <- function(hpdTrees) {
  data = hpdTrees %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    group_by(tree) %>% 
    summarise(prop = sum(recovery) * 100 / nSim)
  
  g = ggplot(data, aes(x = tree, y = prop, color = tree)) +
    geom_point(size = 2) +
    geom_linerange(aes(ymin = 0, ymax = prop), linewidth = 0.1) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    scale_x_discrete(labels = NULL) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Tree coverage (%)", color = NULL) 
  
  return(g)
}


plot_size_boxplot <- function(hpdTrees) {
  data = hpdTrees %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees),
           size = nTrees / nTotal) %>%
    select(tree, seed, size)
  
  g = ggplot(data, aes(x = tree, y = size, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_discrete(labels = NULL) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative size of\n95% credible tree set", color = NULL) 
  
  return(g)
}


# alternative metrics for MCC trees
plot_shPI_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, shPI)
  
  g = ggplot(data, aes(x = tree, y = shPI, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) + 
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared PI", color = NULL)
  
  return(g)
}


plot_Nye_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, Nye)
  
  g = ggplot(data, aes(x = tree, y = Nye, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) + 
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Nye similarity", color = NULL)
  
  return(g)
}


plot_KS_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, KS)
  
  g = ggplot(data, aes(x = tree, y = KS, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) + 
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "KS distance", color = NULL)
  
  return(g)
}


plot_WS_boxplot <- function(mccP) {
  data = mccP %>% 
    bind_rows(.id = 'tree') %>%
    mutate(tree = factor(tree, levels = trees)) %>%
    select(tree, seed, WS)
  
  g = ggplot(data, aes(x = tree, y = WS, color = tree)) +
    geom_boxplot(show.legend = FALSE, outlier.size = 0.5) + 
    scale_x_discrete(labels = NULL) +
    scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2, 0.5)) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Wasserstein distance", color = NULL)
  
  return(g)
}


plot_summary <- function(path, title) {
  hpdTrees = readRDS(file.path(path, "hpdTreesOutput.Rdat"))
  mccP = readRDS(file.path(path, "mccOutput.Rdat"))
  
  g1 = plot_shPI_boxplot(mccP)
  g2 = plot_Nye_boxplot(mccP)
  g3 = plot_KS_boxplot(mccP)
  g4 = plot_WS_boxplot(mccP)
  g5 = plot_coverage_lollipop(hpdTrees)
  g6 = plot_size_boxplot(hpdTrees)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3 + g4 + g5 + g6 +
    plot_layout(nrow = 1)
  return(g)
}

g_ns = plot_summary(paths['non-sequential'], 'A: non-sequential recordings') 
g_s = plot_summary(paths['sequential'], 'B: sequential recordings')

g_ns / g_s +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")
ggsave("pdf/SuppFigure3_BaselineTreeInference.pdf", height = 8, width = 11)

