# -----
# Script for filtering out missing barcodes from simulated alignments
# -----

library(stringr)
library(dplyr)
library(ape)

codeDir = "~/Projects/celldev"
dataDir = "~/Projects/celldev_data"

# set path to repo as working directory
setwd(codeDir)

# load functions for reading alignments/ tapes
source("analysis_functions.R")

# SETUP
nSim = 20
trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")
rho = 0.1 # sampling proportion in case of sampling
method = "TiDe" # TiDe or Typewriter
nBarcodes = 20
rawDir = paste0(dataDir, "/homog/", method, "/baseline_noise/simulationOutputRaw")
simDir = paste0(dataDir, "/homog/", method, "/baseline_noise/simulationOutput")
if (!dir.exists(simDir)) { dir.create(simDir)}

# filtering threshold: the maximal number of targets with at least this proportion of retained cells will be selected
threshold = 0.1

# all simulations
runs = with(expand.grid(tree = trees, sim = seq_len(nSim)), paste0(tree, "_", sim))

# get raw alignments
if (method == "TiDe") {
  alignments = sapply(runs, function(x) {
    f = paste0(rawDir, "/", x, ".alignment.nexus") # get raw alignment file
    data = load_alignment(f) 
    data[data == 51] = NA # assign missing label
    return(data)
  }, simplify = F, USE.NAMES = T)
  
} else if (method == "Typewriter") {
  alignments = sapply(runs, function(x) {
    # get all tapes
    fs = paste0(rawDir, "/", x, ".alignment_", c(1:nBarcodes), ".nexus")
    tapes = load_tapes(fs, num = F)
    
    # summarized alignment (either number of edits per tape or NA)
    data = lapply(c(1:nBarcodes), function(i) {
      tapeID = paste0("T", i, "_")
      tape = tapes %>% select(contains(tapeID))
      missing = sapply(c(1:nrow(tape)), function(c) { 
        if (all(tape[c, ] == "?")) {
          return(NA)
        } else {
          return(sum(as.numeric(tape[c, ]) > 0))
        } 
      })
      df = data.frame(missing)
      rownames(df) = rownames(tape)
      colnames(df) = paste0("T", i, "")
      return(df)
    }) %>% bind_cols()
    
    return(data)
  }, simplify = F, USE.NAMES = T)
}

  
# look at proportion of missing data (~50%)
sapply(alignments, function(alignment) { sum(is.na(alignment)) / (nrow(alignment)*ncol(alignment)) }) %>% summary


# filter alignments to only include cells with non-missing targets
filteredAlignments = lapply(alignments, function(alignment) {
  
  # total number of cells
  n_cells = nrow(alignment)
  
  # sort targets alignment# sort targets by number of non-missing cells 
  cells_per_target = colSums(!is.na(alignment)) %>% sort(decreasing = T)
  targets = names(cells_per_target)
  
  # iteratively retain targets 
  retained_cells = numeric(length(targets))
  subsets = vector("list", length(targets))
  for (k in seq_along(targets)) {
    current_targets = targets[1:k]
    complete_idx = rowSums(!is.na(alignment[, current_targets, drop = FALSE])) == k
    retained_cells[k] = sum(complete_idx)
    subsets[[k]] = alignment[complete_idx, current_targets, drop = FALSE]
  }
  
  res = data.frame(num_targets = seq_along(targets), retained_cells = retained_cells) %>%
    mutate(retained_cells_prop = retained_cells / n_cells)
  
  # find maximal number of targets with sufficient cells
  res_max = res %>% filter(retained_cells >= n_cells * threshold) %>% slice_tail(n = 1)
  subset_max = subsets[[res_max$num_targets]]
  # add initial tax label for reading .nexus alignment in .xml
  res_max$tax_label = rownames(subset_max)[1]
  
  return(list(stats = res_max, data = subset_max))
})


# collect statistics 
filteredStats = lapply(filteredAlignments, "[[", 1) %>% bind_rows(.id = "tree")
summary(filteredStats$num_targets)
summary(filteredStats$retained_cells)
summary(filteredStats$retained_cells_prop)

# calculate resulting sampling proportion
filteredStats = filteredStats %>% 
  mutate(tree_model = str_extract(tree, "tree_[a-z]+")) %>%
  mutate(seed = str_extract(tree, "[0-9]+")) %>%
  mutate(sampling_prop = ifelse(tree_model %in% c("tree_ss", "tree_sds", "tree_bd"), 
                                retained_cells_prop * rho, retained_cells_prop)) # multiply by rho


# save table for inference
write.csv(filteredStats, file = file.path(simDir, "filtering_stats.csv"), row.names = F, quote = F)


# write filtered alignments to .nexus files
all(names(filteredAlignments) == runs)

if (method == "TiDe") {
  # simply write all alignments
  lapply(runs, function(x) {
    alignment = filteredAlignments[[x]]$data
    file = paste0(simDir, "/", x, "alignment.nexus")
    write.nexus(alignment, file)
  })
  
} else if (method == "Typewriter") {
  lapply(runs, function(x) {
    subset = filteredAlignments[[x]]$data
    valid_cells = rownames(subset)
    valid_tapes = str_extract(colnames(subset), "[0-9]+")
    i = 1
    for (tape in valid_tapes) {
      f_in = paste0(rawDir, "/", x, ".alignment_", tape, ".nexus") # get raw alignment file
      data = load_alignment(f_in, num = F) 
      data = data[valid_cells, ]
      all(data != "?") # check missingness
      f_out = paste0(simDir, "/", x, ".alignment_", i, ".nexus")
      write.nexus(data, f_out)
      i = i + 1
    }
  })
}



# FUNCTIONS
write.nexus <- function(alignment, file) {
  ntax = nrow(alignment)
  taxlabels = rownames(alignment)
  nchar = ncol(alignment)
  alignment_lines = sapply(c(1:ntax), function(i) {
    line = paste0(taxlabels[i], " ", paste(alignment[i, ], collapse = ","))
    if (i == ntax) { line = paste0(line, ";") }
    return(line)
  })
  
  # header
  cat(sprintf(r"(#NEXUS

begin taxa;
	dimensions ntax=%d;
	taxlabels %s;
end;

begin characters;
	dimensions nchar=%d;
	format datatype=editData;
	matrix
)",
ntax,
paste(taxlabels, collapse = " "),
nchar
  ), file = file)
  
  # matrix
  cat(paste0("\t\t", alignment_lines, "\n"), sep = "", file = file, append = T)
  cat("end;\n", file = file, append = T)
}


### Important: also store filtered trees and their parameters!
treeDirIn = paste0(dataDir, "/homog/", "LargeTrees")
treeDirOut = paste0(dataDir, "/homog/", method, "/baseline_noise/filteredTrees")
if (!dir.exists(treeDirOut)) { dir.create(treeDirOut)}

trees = sapply(runs, function(x) {
  f = paste0(treeDirIn, "/", x, ".newick") 
  tree = read.tree(f)
  tree = keep.tip(tree, filteredAlignments[[x]]$data %>% rownames)
  return(tree)
}, simplify = F, USE.NAMES = T)
# write filtered trees
lapply(runs, function(x) { write.tree(trees[[x]], file = paste0(treeDirOut, "/", x, ".newick")) })

treeParams = read.csv(paste0(treeDirIn, "/large_tree_params.csv"))
# align order
trees = trees[treeParams$tree]
# update true tree parameters
treeParams$nCells = sapply(trees, Ntip)
treeParams$rho = filteredStats[match(treeParams$tree, filteredStats$tree), "sampling_prop"]
treeParams$treeLength = sapply(trees, function(t) sum(t$edge.length))
treeParams$treeHeight = sapply(trees, function(t) max(node.depth.edgelength(t)))
# save updated table
write.csv(treeParams, paste0(treeDirOut, "/filtered_tree_params.csv"), row.names = F, quote = F)
