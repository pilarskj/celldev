# -----
# Script for finding a good sequential editing rate
# -----

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(assertthat)

# set path to repo as working directory
# setwd("~/Projects/celldev")
source('analysis_functions.R')

data_dir = '~/Projects/celldev_data'

# simulated Typewriter 1x20 barcodes with different editing rate (0.5, 0.4, 0.45,...)
# comparing the sequential to non-sequential (TiDe baseline) alignments:
td_dir = paste0(data_dir, '/TiDe/baseline/simulationOutput')
tw_dir = paste0(data_dir, '/Typewriter/sequential/simulationOutput')

metrics = c("nEdited", "nEdits", "nFull", "nBarcodes", "nCells", "barcodeDiv", "distHamming")

# summarize statistics on the simulated alignments
simStats = lapply(c(td_dir, tw_dir), function(d) {
  
  # get alignments
  alignments = sapply(list.files(d, full.names = T), load_alignment, simplify = F) 
  names(alignments) = names(alignments) %>% str_extract("tree_.*_[0-9]+") %>% sub("\\..*", "", .)
  
  # get unique barcodes per alignment
  barcodes = lapply(alignments, function(alignment) {
    alignment %>% group_by_all() %>% summarize(count=n()) %>% arrange(desc(count)) %>% 
      suppressMessages() })
  
  # set up stats
  stats = data.frame(matrix(nrow = length(alignments), ncol = 8))
  colnames(stats) = c("tree", metrics)
  stats$tree = names(alignments)
  
  # get median number of edited targets/tapes per cell 
  stats$nEdited = sapply(alignments, function(a) {median(rowSums(a != 0)) })

  # get number of different edits made across targets and cells
  stats$nEdits = sapply(alignments, function(a) {length(unique(unlist(a)[unlist(a) != 0])) })
  
  # get number of 'full' targets
  stats$nFull = sapply(alignments, function(a) {sum(rowSums(a == 20))})
  
  # get number of cells
  stats$nCells = sapply(alignments, nrow)
  
  # get number of unique barcodes 
  stats$nBarcodes = sapply(barcodes, nrow)
  
  # calculate barcode diversity: number of barcodes per cell (how many cells share a barcode?)
  stats$barcodeDiv = stats$nBarcodes / stats$nCells 
  
  # calculate median pairwise Hamming distance
  stats$distHamming = sapply(alignments, calc_hamming_dist)
  
  return(stats)
}) 
names(simStats) = c('TiDe_baseline', 'Typewriter_sequential')


# check and reformat tables for plotting
assert_that(all(simStats$TiDe_baseline$nCells == simStats$Typewriter_sequential$nCells))
compStats = simStats %>%
  bind_rows(.id = 'setting') %>%
  pivot_longer(cols = all_of(metrics), names_to = 'metric', values_to = 'value') %>%
  filter(metric != 'nCells')


# plot comparison of the barcode metrics
ggplot(compStats, aes(x = setting, y = value)) +
  geom_boxplot() +
  expand_limits(y = 0) +
  labs(x = NULL, y = NULL) +
  facet_wrap(vars(metric), scales = 'free_y')

# TODO: apply statistical test?
