# Results for noisy data

setwd("~/Projects/celldev/figures")
source("homog_functions.R")

# also report growth rate
params = c('birthRate', 'deathRate', 'growthRate', 'editRate', 'treeHeight', 'treeLength') 
params_labels = c(params_labels, 'growthRate' = 'growth rate')


collect_plots <- function(method, title) {
  
  if (method == "TiDe") {
    scenarios = c("nTargets_5", "baseline_noise")
    sc_labels = c("nTargets_5" = "5 targets", "baseline_noise" = "filtered\ntargets")
  } else if (method == "Typewriter") {
    scenarios = c("nTapes_5", "baseline_noise")
    sc_labels = c("nTapes_5" = "5 tapes", "baseline_noise" = "filtered\ntapes")
  }
  
  paths = sapply(scenarios, function(x) paste0(dir, "/", method, "/", x, "/analysisOutput"))
  
  # for scenario with noise, show coverage and bias
  infP = readRDS(file.path(paths['baseline_noise'], "/infPerformance.Rdat"))
  
  cov = lapply(infP, '[[', 'coverage') %>%
    bind_rows(.id = 'tree') %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'coverage') %>%
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = params)) %>%
    group_by(tree, parameter) %>%
    summarise(prop = sum(coverage) * 100 / nSim)
   
  # gCov = ggplot(cov, aes(x = parameter, y = prop, color = tree)) +
  #   geom_point(size = 2, position = position_dodge(width = 0.6)) + 
  #   geom_linerange(aes(ymin = 0, ymax = prop), position = position_dodge(width = 0.6), linewidth = 0.1) +
  #   geom_hline(yintercept = 80, linetype = 'dashed') + # threshold at 80%
  #   scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  #   scale_x_discrete(labels = params_labels) +
  #   scale_color_manual(values = palette, labels = tree_labels) +
  #   labs(x = NULL, y = "Coverage (%)", color = NULL) +
  #   theme(legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")
  # 
  # legend <<- get_legend(gCov)
  
  bias = lapply(infP, '[[', 'bias_rel') %>% 
    bind_rows(.id = 'tree') %>%
    pivot_longer(cols = all_of(params), names_to = 'parameter', values_to = 'bias') %>%
    mutate(tree = factor(tree, levels = trees), parameter = factor(parameter, levels = params))
  
  gBias = ggplot(bias, aes(x = parameter, y = bias, color = tree)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_hline(yintercept = 0, linetype = 'dashed') + # line at 0
    ylim(-1, 2) +
    scale_x_discrete(labels = params_labels) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Relative bias", color = NULL) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
          legend.position = 'bottom', legend.direction = "vertical", legend.justification = "right")
  
  legend <<- get_legend(gBias)
  
  # tree inference (compared to the limited capacity scenario)
  mccP = sapply(paths, function(p) readRDS(file.path(p, 'mccOutput.Rdat')), simplify = F, USE.NAMES = T)
  mcc = lapply(mccP, function(x) {x %>% 
      bind_rows(.id = 'tree') }) %>%
    bind_rows(.id = 'scenario') %>%
    mutate(tree = factor(tree, levels = trees),
           scenario = factor(scenario, levels = scenarios))
  
  g1 = ggplot(mcc, aes(x = scenario, y = wRF, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    scale_x_discrete(labels = sc_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Weighted RF distance", color = NULL)
  
  g2 = ggplot(mcc, aes(x = scenario, y = shPI, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    scale_x_discrete(labels = sc_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "Shared PI", color = NULL)
  
  g3 = ggplot(mcc, aes(x = scenario, y = KS, color = tree)) + 
    geom_boxplot(outlier.size = 0.5, show.legend = F) + 
    scale_x_discrete(labels = sc_labels) +
    scale_y_continuous(limits = c(0, 1), breaks = scales::pretty_breaks()) +
    scale_color_manual(values = palette, labels = tree_labels) +
    labs(x = NULL, y = "KS distance", color = NULL)
  
  g = gBias + ggtitle(title) + theme(plot.title.position = "plot", legend.position = "none") + g1 + g2 + g3 +
    plot_layout(widths = c(1.5, 1, 1, 1)) &
    theme(axis.text.y = element_text(size = 8))

  return(g)
}


g_ns = collect_plots("TiDe", "A: non-sequential recordings (noisy)")
g_s = collect_plots("Typewriter", "B: sequential recordings (noisy)")

(wrap_elements(g_ns) / wrap_elements(g_s)) / wrap_elements(legend) +
  plot_layout(heights = c(1, 1, 0.25)) 
ggsave("pdf/Figure5_Noise.pdf", height = 8, width = 10)

