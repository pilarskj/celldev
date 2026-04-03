# Figure with type statistics

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(showtext)
library(patchwork)
library(HDInterval)

# plotting options
palette = c("0" = "#56B4E9", "1" = "#009E73", "2" = "#E69F00", "3" = "#CC79A7")
font_add("lmroman", regular = "~/lmroman10-regular.otf") # Latex font
showtext_auto()
theme_set(theme_classic(base_size = 12, base_family = 'lmroman') +
            theme(plot.title.position = "plot"))

# paths
data_dir = "/Volumes/stadler/cEvoUnpublished/2023-Julia-Celldev/heterog"
code_dir = "~/Projects/celldev"
setwd(file.path(code_dir, "figures"))

# for reversing time
Texp = 40

plot_typestats <- function(method, transition, subtitle) {
  
  # true stats
  trueStats = read.csv(paste0(data_dir, "/Trees/", transition, "/transition_stats.csv"))
  seeds = sort(unique(trueStats$tree)) %>% as.character
  ids = c(1:length(seeds)) %>% setNames(seeds)
  trueStats = trueStats %>%
    mutate(tree = as.character(tree)) %>%
    mutate(ID = ids[tree])
  
  # posterior stats
  postStats = readRDS(paste0(data_dir, "/", method, "/", transition, "/analysisOutput/treeStats.Rdat"))
  
  # combine
  data = postStats %>%
    left_join(trueStats %>%
                rename_with(.fn = ~paste0(.x, "_true"), .cols = c(time_spent, n_trans, first_trans)),
              by = c("tree", "type")) %>%
    pivot_longer(
      cols = -c(tree, ID, type, n_tips), 
      names_to = c("parameter", ".value"), 
      names_pattern = "(.*)_(lower|median|upper|true)"
    ) %>%
    mutate(covered = ifelse(true >= lower & true <= upper, T, F)) %>%
    mutate(type_cov = ifelse(covered, type, NA) %>% as.factor) # for plotting
  
  g1 = ggplot(data %>% filter(type != 0, parameter == "n_trans"), aes(x = ID)) +
    geom_pointrange(aes(ymin = lower, y = median, ymax = upper, color = type_cov), size = 0.1, alpha = 0.6, linewidth = 1, show.legend = F) + 
    geom_point(aes(y = true), size = 0.7) +
    expand_limits(y = 0) +
    scale_color_manual(values = palette, na.value = "grey") +
    facet_wrap(~type, labeller = as_labeller(function(x) paste("to type", x))) +
    labs(x = "Simulation", y = "Number of transitions") 
  
  g2 = ggplot(data %>% filter(type != 0, parameter == "first_trans"), aes(x = ID)) +
    geom_pointrange(aes(ymin = Texp - lower, y = Texp - median, ymax = Texp - upper, color = type_cov), size = 0.1, alpha = 0.6, linewidth = 1, show.legend = F) + 
    geom_point(aes(y = Texp - true), size = 0.7) +
    ylim(c(0, 40)) + 
    scale_color_manual(values = palette, na.value = "grey") +
    facet_wrap(~type, labeller = as_labeller(function(x) paste("to type", x))) +
    labs(x = "Simulation", y = "Time of first transition") 
  
  g3 = ggplot(data %>% filter(parameter == "time_spent"), aes(x = ID)) +
    geom_pointrange(aes(ymin = lower, y = median, ymax = upper, color = type_cov), size = 0.1, alpha = 0.6, linewidth = 1, show.legend = F) + 
    geom_point(aes(y = true), size = 0.7) +
    expand_limits(y = 0) +
    scale_color_manual(values = palette, na.value = "grey") +
    facet_wrap(~type, nrow = 1, labeller = as_labeller(function(x) paste("in type", x))) +
    labs(x = "Simulation", y = "Sum of branch lengths") 
  
  g = (g1 + labs(subtitle = subtitle)) / g2 / g3 + plot_layout(axis_titles = "collect") &
    theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
  
  return(g)
}

p1 = plot_typestats("TiDe", "distinct", "(I) terminal transitions")
p2 = plot_typestats("TiDe", "hierarchical", "(II) chain-like transitions")
p3 = plot_typestats("Typewriter", "distinct", "(I) terminal transitions")
p4 = plot_typestats("Typewriter", "hierarchical", "(II) chain-like transitions")

p_ns = (p1 | p2) + plot_layout(axis_titles = "collect") + plot_annotation(title = "A: non-sequential recordings")
p_s = (p3 | p4) + plot_layout(axis_titles = "collect") + plot_annotation(title = "B: sequential recordings")

pdf('pdf/SuppFigure16_TypeStats.pdf', height = 10, width = 9) 
p_ns
p_s
dev.off()

