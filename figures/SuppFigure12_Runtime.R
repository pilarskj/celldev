# Graphs of runtime (recorded for new scenarios)

library(tracerer)

setwd("~/Projects/celldev/figures")
source("/homog_functions.R") # read in global settings

settings = c("baseline_large", "baseline_sampled", "baseline_noise")


# function for converting times
convert_to_minutes <- function(x) {
  h <- as.numeric(str_extract(x, "\\d+(?=h)"))
  m <- as.numeric(str_extract(x, "\\d+(?=m)"))
  s <- as.numeric(str_extract(x, "\\d+(?=s)"))
  
  h[is.na(h)] <- 0
  m[is.na(m)] <- 0
  s[is.na(s)] <- 0
  
  return(60*h + m + s/60)
}


collect_runtimes <- function(method, settings) {
  
  # read recorded runtimes
  paths = sapply(settings, function(x) paste0(dir, "/", method, "/", x, "/inferenceOutput/runtimes.csv"))
  df = lapply(paths, read.csv) %>% bind_rows(.id = "setting")
  
  # convergence (posterior ESS > 200)
  df$chainLength = NA
  df$nSamples = NA
  df$converged = NA
  for (i in c(1:nrow(df))) {
    setting = df[i, "setting"]
    tree = df[i, "tree"]
    
    # for each setting and for each tree, read log file, remove 10% burn-in
    file = paste0(dir, "/", method, "/", setting, "/inferenceOutput/", tree, ".inference.log")
    log = remove_burn_ins(parse_beast_tracelog_file(file), burn_in_fraction = 0.1)
    
    # get final chain length
    df[i, "chainLength"] = log[nrow(log), "Sample"]
    
    # estimate how ESS grows wrt. samples collected
    ess = lapply(seq(100, nrow(log), 100), function(k) {
      data.frame(nSamples = log[k, "Sample"], ESS = calc_ess(log[1:k, "posterior"], sample_interval = 5000)) }) %>% 
      bind_rows()
    
    # store number of samples needed for convergence
    if (ess[nrow(ess), "ESS"] > 200) {
      df[i, "nSamples"] = ess %>% filter(ESS > 200) %>% slice_head(n = 1) %>% pull(nSamples)
      df[i, "converged"] = TRUE
    } else {
      # for runs which did not converge, extrapolate chainLength
      finalESS = calc_ess(log$posterior, sample_interval = 5000)
      df[i, "nSamples"] = chainLength / finalESS * 200
      df[i, "converged"] = FALSE
    }
    if (i %% 10 == 0) { # for progress tracking
      print(paste0("row ", i, " done"))
    }
  }
  
  # convert times
  df = df %>% 
    mutate(runtime_min = convert_to_minutes(runtime)) %>% # Msamples/min
    mutate(runtime_chain = chainLength / 1e6 * runtime_min / 60, # total/h
           runtime_conv = nSamples / 1e6 * runtime_min / 60) # convergence/h
  
  # add tree sizes
  tree_paths = c("baseline_large" = paste0(dir, "/LargeTrees/tree_params.csv"),
                 "baseline_sampled" = paste0(dir, "/SampledTrees/tree_params.csv"),
                 "baseline_noise" = paste0(dir, "/", method, "/baseline_noise/filteredTrees/filtered_tree_params.csv"))
  tree_params = lapply(tree_paths, read.csv) %>% bind_rows(.id = "setting")
  df = df %>% left_join(tree_params %>% select(setting, tree, nCells), by = c("tree", "setting"))
  
  # format columns
  df_tw = df_tw %>% 
    mutate(seed = as.numeric(str_extract(tree, '[0-9]+')),
           tree = factor(str_extract(tree, 'tree_[a-z]+'), levels = trees),
           method = method) %>%
    arrange(setting, tree, seed)
  
  return(df)
}


# run this (slow)
# df_td = collect_runtimes("TiDe", settings)
# df_tw = collect_runtimes("Typewriter", settings)
# df = rbind(df_td, df_tw)
# write.csv(df, "runtimes.csv", quote = F, row.names = F)
# df = read.csv("runtimes.csv") %>% mutate(tree = factor(tree, trees))

# now plot

method_labels = c("TiDe" = "TiDeTree", "Typewriter" = "SciPhy")

# runtime wrt. tree size
g1 = ggplot(df, aes(x = nCells, y = runtime_conv, color = tree, alpha = converged)) +
  geom_point(size = 0.8) +
  scale_x_continuous(limits = c(0, 2000), breaks = pretty_breaks()) +
  expand_limits(y = 100) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_manual(values = palette, labels = tree_labels) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.25)) +
  labs(x = "Tree size (#sampled cells)", y = "Time until convergence (h)", color = NULL) +
  facet_wrap(~method, labeller = as_labeller(method_labels), scales = "free")

# runtime of method grouped by tree size
df = df %>% 
  mutate(treeSize = ifelse(setting == "baseline_large", "large", "small")) %>%
  mutate(treeSize = factor(treeSize, levels = c("small", "large")))
# large: 332-1668 cells, small: 33-250

size_labels = c("small" = "small trees (<250 tips)",
                "large" = "large trees (>250 tips)")

# time per Msamples
factor_min_df = df %>% select(method, setting, tree, seed, runtime = runtime_min, treeSize) %>%
  group_by(treeSize) %>% mutate(max_runtime = quantile(runtime, probs = 0.995)) %>%
  pivot_wider(names_from = "method", values_from = "runtime") %>%
  mutate(factor = TiDe / Typewriter) %>%
  group_by(treeSize) %>%
  summarise(factor = median(factor), max_runtime = unique(max_runtime)) %>%
  mutate(group1 = "TiDe", group2 = "Typewriter", y.position = max_runtime, label = paste0("x", round(factor, 0)))

g2 = ggplot(df, aes(x = method, y = runtime_min, color = tree)) +
  geom_boxplot(outlier.size = 0.5, show.legend = F) +
  stat_pvalue_manual(data = factor_min_df) +
  scale_x_discrete(labels = method_labels) +
  expand_limits(y = 0) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "Runtime per 1M samples (min)", color = NULL) +
  facet_wrap(~treeSize, labeller = as_labeller(size_labels), scales = "free_y")

factor_conv_df = df %>% select(method, setting, tree, seed, runtime = runtime_conv, treeSize) %>%
  group_by(treeSize) %>% mutate(max_runtime = quantile(runtime, probs = 0.995)) %>%
  pivot_wider(names_from = "method", values_from = "runtime") %>%
  mutate(factor = TiDe / Typewriter) %>%
  group_by(treeSize) %>%
  summarise(factor = median(factor), max_runtime = unique(max_runtime)) %>%
  mutate(group1 = "TiDe", group2 = "Typewriter", y.position = max_runtime, label = paste0("x", round(factor, 0)))

g3 = ggplot(df, aes(x = method, y = runtime_conv, color = tree)) +
  geom_boxplot(outlier.size = 0.5, show.legend = F) +
  stat_pvalue_manual(data = factor_conv_df) +
  scale_x_discrete(labels = method_labels) +
  expand_limits(y = 0) +
  scale_color_manual(values = palette) +
  labs(x = NULL, y = "Time until convergence (h)", color = NULL) +
  facet_wrap(~treeSize, labeller = as_labeller(size_labels), scales = "free_y")

# combine plots
(g1 + ggtitle("A")) / (g2 + ggtitle("B")) / g3 & theme(plot.title.position = "plot")
ggsave("pdf/SuppFigure12_Runtime.pdf", height = 10, width = 10)
