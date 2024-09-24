library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

setwd("~/Desktop/Master Thesis/Code and Data/Part1")
source('analysis_functions.R')

tdt_dir = '~/Desktop/Master Thesis/Code and Data/Part1/TiDe/baseline/simulationOutput'
tw_dir = '~/Desktop/Master Thesis/Code and Data/Part1/Typewriter/sequential/simulationOutput'

# load alignments
tdt_sim = sapply(list.files(tdt_dir, full.names = T), load_alignment, simplify = F) 
names(tdt_sim) = str_extract(names(tdt_sim), "tree_.*_[0-9]+")

tw_sim = sapply(list.files(tw_dir, full.names = T), load_alignment, simplify = F) 
names(tw_sim) = str_extract(names(tw_sim), "tree_.*_[0-9]+")


# ALIGNMENTS STATS
# get median number of edited targets/tapes per cell 
summary(sapply(tdt_sim, function(a) {median(rowSums(a != 0))}))
summary(sapply(tw_sim, function(a) {median(rowSums(a != 0))}))

# get number of different edits made across targets and cells
summary(sapply(tdt_sim, function(a) {length(unique(unlist(a)[unlist(a) != 0])) }))
summary(sapply(tw_sim, function(a) {length(unique(unlist(a)[unlist(a) != 0])) }))

# get number of 'full' barcodes
summary(sapply(tdt_sim, function(a) {sum(rowSums(a == 20))}))
summary(sapply(tw_sim, function(a) {sum(rowSums(a == 20))}))


# BARCODE STATS
# get unique barcodes per alignment
tdt_bc = lapply(tdt_sim, function(alignment) {
  alignment %>% group_by_all() %>% summarize(count=n()) %>% arrange(desc(count)) %>% 
    suppressMessages() })

tw_bc = lapply(tw_sim, function(alignment) {
  alignment %>% group_by_all() %>% summarize(count=n()) %>% arrange(desc(count)) %>% 
    suppressMessages() })

# get number of unique barcodes 
summary(sapply(tdt_bc, nrow))
summary(sapply(tw_bc, nrow))

# calculate barcode diversity: number of barcodes per cell (how many cells share a barcode?)
summary(sapply(tdt_bc, nrow) / sapply(tdt_sim, nrow))
summary(sapply(tw_bc, nrow) / sapply(tw_sim, nrow))




# Hamming distances
calc_hamming_dist <- function(x) { #alignment) {
  # hdist = matrix(nrow = nrow(alignment), ncol = nrow(alignment))
  # for (i in c(1:(nrow(alignment)-1))) {
  #   for (j in c((i+1):nrow(alignment))) {
  #     hdist[i, j] = sum(alignment[i, ] != alignment [j, ])
  #   }
  # }
  xt <- t(x)
  hdist = sapply(1:nrow(x), function(y) colSums(xt != xt[, y]))
  return(median(hdist, na.rm = T))
}

summary(sapply(tdt_bc, calc_hamming_dist))
summary(sapply(tw_bc, calc_hamming_dist))

