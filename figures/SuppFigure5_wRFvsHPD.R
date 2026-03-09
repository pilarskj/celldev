# Association between wRF, branch uncertainty, and relative HPD width of phylodynamic parameters

library(treeio)

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

scenarios = c("non-sequential", "sequential")
paths = c("non-sequential" = paste0(dir, "/TiDe/baseline"), "sequential" = paste0(dir, "/Typewriter/baseline"))

runs = with(expand.grid(tree = trees, sim = seq_len(nSim)), paste0(tree, "_", sim))

plot_association <- function(path, title) {
  infP = readRDS(paste0(path, "/analysisOutput/infPerformance.Rdat"))
  mccP = readRDS(paste0(path, "/analysisOutput/mccOutput.Rdat"))
  mccTrees = sapply(runs, function(x) {
    f = paste0(path, "/inferenceOutput/", x, ".mcc.tree")
    tree = read.beast(f)
    return(tree)
  }, simplify = F, USE.NAMES = T)
  
  # extract relative HPD width
  hpd = lapply(infP, '[[', 'hpd_width_rel') %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, all_of(params))

  # extract wRF distance
  wRF = mccP %>% 
    bind_rows(.id = 'tree') %>%
    select(tree, seed, wRF) %>%
    mutate(seed = as.character(seed))
  
  # evaluate branch uncertainty in summary tree
  branchUnc = lapply(mccTrees, function(mcc) {
    df = as_tibble(mcc) 
    df = df %>% 
      filter(is.na(label)) %>% # only internal nodes
      select(node, height_median, height_0.95_HPD) %>%
      filter(!is.na(height_median)) %>%
      mutate(height_lower = map_dbl(height_0.95_HPD, 1), height_upper = map_dbl(height_0.95_HPD, 2)) %>%
      mutate(branch_unc = (height_upper - height_lower) / height_median) # relative HPD width around internal nodes
    return(df)
  }) %>% bind_rows(.id = 'tree_id') %>%
    group_by(tree_id) %>%
    summarise(mean_branch_unc = mean(branch_unc), sd_branch_unc = sd(branch_unc)) %>%
    mutate(tree = str_extract(tree_id, "tree_[a-z]+"), seed = str_extract(tree_id, "[0-9]+")) %>%
    select(-tree_id)
  
  data = hpd %>%
    full_join(wRF, by = c("tree", "seed")) %>%
    full_join(branchUnc, by = c("tree", "seed")) %>%
    mutate(tree = factor(tree, levels = trees))
  
  # g1 = ggplot(data, aes(x = wRF, y = birthRate)) +
  #   geom_point(aes(color = tree)) +
  #   stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.1, color = "black") +
  #   stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
  #   expand_limits(x = c(0,0.2), y = c(0,1)) +   
  #   scale_y_continuous(breaks = pretty_breaks()) +
  #   scale_color_manual(values = palette, labels = tree_labels) +
  #   labs(x = "Weighted RF distance", y = "Relative HPD width\n(cell division rate)", color = NULL)
  
  g1 = ggplot(data, aes(x = wRF, y = mean_branch_unc)) +
    geom_pointrange(aes(color = tree, ymin = mean_branch_unc - sd_branch_unc/2, ymax = mean_branch_unc + sd_branch_unc/2),
                    size = 0.1, linewidth = 0.05) +
    stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.05, color = "black") +
    stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
    expand_limits(x = 0, y = 0.5) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = "Weighted RF distance", y = "Relative HPD width\n(branching times)", color = NULL)
  
  g2 = ggplot(data, aes(x = growthRate, y = mean_branch_unc)) +
    geom_pointrange(aes(color = tree, ymin = mean_branch_unc - sd_branch_unc/2, ymax = mean_branch_unc + sd_branch_unc/2),
                    size = 0.1, linewidth = 0.05) +
    stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.05, color = "black") +
    stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
    expand_limits(x = c(0,1), y = 0.5) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = "Relative HPD width\n(growth rate)", y = "Relative HPD width\n(branching times)", color = NULL)
  
  g3 = ggplot(data, aes(x = birthRate, y = mean_branch_unc)) +
    geom_pointrange(aes(color = tree, ymin = mean_branch_unc - sd_branch_unc/2, ymax = mean_branch_unc + sd_branch_unc/2),
                    size = 0.1, linewidth = 0.05) +
    stat_smooth(method = "lm", linewidth = 0.2, alpha = 0.05, color = "black") +
    stat_cor(aes(label = ..r.label..), method = "kendall", cor.coef.name = "cor.") +
    expand_limits(x = c(0,1), y = 0.5) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = "Relative HPD width\n(cell division rate)", y = "Relative HPD width\n(branching times)", color = NULL)
  
  g = g1 + ggtitle(title) + theme(plot.title.position = "plot") +
    g2 + g3

  return(g)
}

# plot figures
g_ns = plot_association(paths['non-sequential'], 'A: non-sequential recordings') 
g_s = plot_association(paths['sequential'], 'B: sequential recordings')

g_ns / g_s +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")
ggsave('pdf/SuppFigure5_wRFvsHPD.pdf', height = 9, width = 10)
