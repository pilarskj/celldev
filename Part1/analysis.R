# ----------
# Script for the analysis of simulated alignments from .nexus files and
# BEAST2 inference results from .log, .hpd.trees and .mcc.tree files

# During inference, the following parameters were estimated from the simulated data: 
# birth rate, death rate, edit rate, tree length, tree height
# The sampling proportion, origin, editing window and edit probabilities were fixed.
# ----------

# load packages
library(tracerer)
library(HDInterval)
library(ape)
library(treeio)
library(TreeDist) 
library(phangorn)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# clear environment
rm(list = ls()) 
# setwd("~/Desktop/Part1/")

# load helper functions
source("analysis_functions.R")


### SCENARIO
# lineage recorder
method = "Typewriter"

# setup
setup = "baseline"
editRate = 0.05

# number of simulations
nSim = 20

# trees
trees = c("tree_s", "tree_ss", "tree_sd", "tree_sds", "tree_bd")

# tree parameters
treeParams = read.csv("Trees/treeParams.csv")

# inferred parameters
parameters = c("treeHeight", "treeLength", "editRate", 
               "birthRate", "deathRate")#, "deathProb", "growthRate", "turnoverRate")
nParams = length(parameters)

# directory
dir = paste0(method, "/", setup, "/")

# output directory
outDir = paste0(dir, "analysisOutput/")
if (!dir.exists(outDir)) {dir.create(outDir)}



### SIMULATED ALIGNMENTS
# summarize statistics on the simulated alignments
simStats = lapply(trees, function(tree) {
  # get alignment files
  simDir = paste0(dir, "simulationOutput")
  simFiles = list.files(simDir, pattern = paste0(tree, "_"), full.names = T)
  n = length(simFiles)
  
  # load alignments
  if (method == "TiDe"){
    alignments = sapply(simFiles, load_alignment, simplify = F) 
    names(alignments) = str_extract(str_extract(names(alignments), paste0(tree, "_[0-9]+")), "[0-9]+") # extract simulation seed
  } else if (method == "Typewriter"){ 
    alignments = sapply(c(1:nSim), function(i) {
      tape_files = simFiles[str_detect(simFiles, paste0("_", i, "[.]alignment"))]
      data = load_tapes(tape_files)
      return(data)
    }, simplify = F)
    names(alignments) = as.character(c(1:nSim))
  }
  
  # get unique barcodes per alignment
  barcodes = lapply(alignments, function(alignment) {
    alignment %>% group_by_all() %>% summarize(count=n()) %>% arrange(desc(count)) %>% 
      suppressMessages() })
  
  # calculate statistics per alignment
  stats = data.frame(matrix(nrow = nSim, ncol = 6))
  colnames(stats) = c("seed", "nEdited", "nEdits", "nBarcodes", "nCells", "barcodeDiv")
  
  stats$seed = names(alignments)
  
  # get median number of edited targets/tapes per cell 
  if (method == "TiDe"){
    stats$nEdited = sapply(alignments, function(a) {median(rowSums(a != 0)) })
  } else if (method == "Typewriter"){
    stats$nEdited = sapply(alignments, function(a) {median(rowSums(a %>% select(contains("V1")) != 0)) })
  }
  
  # get number of different edits made across targets and cells
  stats$nEdits = sapply(alignments, function(a) {length(unique(unlist(a)[unlist(a) != 0])) })
  
  # get number of unique barcodes 
  stats$nBarcodes = sapply(barcodes, nrow)
  
  # check if number of cells in simulations corresponds to number of tips in tree
  if (method == "TiDe"){
    stats$nCells = sapply(simFiles, get_ntaxa)
  } else if (method == "Typewriter"){
    stats$nCells = sapply(c(1:nSim), function(i) {
      get_ntaxa(simFiles[str_detect(simFiles, paste0("_", i, "[.]alignment_1[.]"))])})
  }
  print(all(stats$nCells == treeParams[treeParams$tree == tree, "nCells"]))
  
  # calculate barcode diversity: number of barcodes per cell (how many cells share a barcode?)
  stats$barcodeDiv = stats$nBarcodes / stats$nCells 
  
  # sort by seed
  stats = stats %>% arrange(as.numeric(seed))
  
  return(stats)
})

# store
names(simStats) = trees
saveRDS(simStats, file = paste0(outDir, "simStats.Rdat"))



### INFERENCE RESULTS
# summarize inference output for each tree
infOutput = lapply(trees, function(tree) {
  # find log files
  infDir = paste0(dir, "inferenceOutput/")
  logFiles = list.files(infDir, pattern = paste0(tree, "_.*.log")) 
  n = length(logFiles)
  
  # create output table
  outputDat = data.frame(matrix(nrow = n,  ncol = 3*nParams + 3)) 
  colnames(outputDat) = c(get_sumstat_names(parameters), "seed", "minESS", "chainLength")
  
  for (i in 1:n){
    # get log file and its simulation seed
    file = logFiles[i]
    outputDat[i, "seed"] = str_extract(file, "[0-9]+")
    
    # parse beast log 
    log = remove_burn_ins(parse_beast_tracelog_file(paste0(infDir, file)), burn_in_fraction = 0.1)
    log = log %>% rename(treeHeight = tree.height, treeLength = tree.treeLength)
    
    # check ESS
    ess = calc_esses(log, sample_interval = 5000)
    outputDat[i, "minESS"] = min(ess)
    
    # get chain length
    outputDat[i, "chainLength"] = log$Sample[nrow(log)]
    
    # calculate missing parameter
    log = log %>% mutate(deathProb = deathRate/ (birthRate + deathRate))
    
    # calculate growth and turnover rate
    log = log %>% mutate(growthRate = birthRate - deathRate,
                         turnoverRate = birthRate/deathRate)
    
    # save summary stats in output data frame
    outputDat[i, 1:(3*nParams)] = get_sumstats_from_log(log, parameters)
    
    # print message
    print(paste0(tree, " file ", i, " done"))
  }
  
  outputDat$seed = as.numeric(outputDat$seed)
  outputDat = outputDat[order(outputDat$seed), ] # order according to simulation seed
  return(outputDat)
})

# store
names(infOutput) = trees
saveRDS(infOutput, file = paste0(outDir, "infOutput.Rdat"))

# infOutput = readRDS(paste0(outDir, "infOutput.Rdat"))


### INFERENCE PERFORMANCE
# filter valid simulations
infOutput = lapply(infOutput, function(x) {x %>% filter(minESS >= 200)}) 
print(lapply(infOutput, nrow))

# extract true simulation parameters
trueParams = treeParams %>% 
  mutate(seed = str_extract(tree, "[0-9]+")) %>%
  mutate(tree = str_remove(tree, "_[0-9]+")) %>%
  mutate(editRate = editRate) %>%
  select(tree, seed, all_of(parameters)) 
trueParams = rbind(trueParams %>% slice(rep(which(tree == "tree_s"), each = 20)) %>% mutate(seed = c(1:20)),
                   trueParams) %>%
  filter(!is.na(seed))

# calculate width of priors
priors = list()
priors$editRate = c(hdi(qlnorm, 0.95, meanlog = -3, sdlog = 1),
                    "median" = qlnorm(0.5, meanlog = -3, sdlog = 1))
priors$birthRate = c(hdi(qlnorm, 0.95, meanlog = -1.5, sdlog = 1),
                     "median" = qlnorm(0.5, meanlog = -1.5, sdlog = 1))
priors$deathRate = c(hdi(qexp, 0.95, rate = 25),
                     "median" = qexp(0.5, rate = 25))
priors_hpd_width = sapply(priors, function(x) {as.numeric(x["upper"] - x["lower"])})


# calculate performance metrics
perf_metrics = c("coverage", "abs_error_rel", "bias_rel", "hpd_width_rel",
                 "hpd_proportion", "abs_error", "bias", "hpd_width")

infPerformance = lapply(trees, function(tree) {
  outputDat = infOutput[[tree]]
  
  # define output table
  perf = sapply(perf_metrics, function(x) {
    df = data.frame(matrix(nrow = nrow(outputDat), ncol = length(parameters), data = NA))
    colnames(df) = c(parameters)
    df = df %>% rownames_to_column(var = "seed")
    return(df) }, simplify = F, USE.NAMES = T)
  
  for (i in 1:nrow(outputDat)) {
    seed = outputDat[i, "seed"]
    
    for(parameter in parameters){
      
      # get true and estimated parameter
      truth = as.numeric(trueParams[trueParams$tree == tree & trueParams$seed == seed, parameter])
      median = outputDat[i, paste0(parameter, "_median")]
      hpd_lower = outputDat[i, paste0(parameter, "_lower")]
      hpd_upper = outputDat[i, paste0(parameter, "_upper")]
      
      # calculate normalization variable
      if (truth != 0) {
        divide = truth
      } else {
        divide = mean(outputDat[, paste0(parameter, "_median")])
      }
      
      # calculate normalized performance metrics
      perf$coverage[seed, parameter] = ifelse(truth >= hpd_lower & truth <= hpd_upper, TRUE, FALSE)
      perf$abs_error_rel[seed, parameter] = abs(median - truth) / divide
      perf$bias_rel[seed, parameter] = (median - truth) / divide
      perf$hpd_width_rel[seed, parameter] = (hpd_upper - hpd_lower) / divide
      
      if (parameter %in% c("editRate", "birthRate", "deathRate")) {
        perf$hpd_proportion[seed, parameter] = (hpd_upper - hpd_lower) / priors_hpd_width[parameter]
      }
      
      # calculate not normalized performance metrics for death rate
      if (parameter == "deathRate") {
        perf$abs_error[seed, parameter] = abs(median - truth)
        perf$bias[seed, parameter] = median - truth
        perf$hpd_width[seed, parameter] = hpd_upper - hpd_lower
      }
    }
  }
  return(perf)
})

# store
names(infPerformance) = trees
saveRDS(infPerformance, file = paste0(outDir, "infPerformance.Rdat"))



### TREE INFERENCE
# load true trees from .newick files
trueTrees = sapply(trees, function(tree) {
  sapply(c(1:nSim), function(i) {
    read.tree(file = paste0("Trees/", tree, "_", i, ".newick"))
  }, simplify = F)
}, simplify = F, USE.NAMES = T)

# analyse 95% credible sets of trees
hpdTreesOutput = lapply(trees, function(tree) {
  # find files with tree HPDs
  infDir = paste0(dir, "inferenceOutput/")
  treesFiles = list.files(infDir, pattern = paste0(tree, "_[0-9]+.hpd.trees")) 
  n = length(treesFiles)
  
  # create output table
  outputDat = data.frame(matrix(nrow = n,  ncol = 5)) 
  colnames(outputDat) = c("seed", "nTrees", "nTotal", "recovery", "pRank1")
  
  for (i in 1:n) {
    # get trees file and its simulation seed
    file = treesFiles[i]
    filePath = paste0(infDir, file)
    
    seed = as.numeric(str_extract(file, "[0-9]+"))
    outputDat[i, "seed"] = seed
    
    trueTree = trueTrees[[tree]][[seed]]
    
    # if file not empty
    if (!is.na(readLines(filePath)[4])) { 
      outputDat[i, "nTotal"] = get_ntrees(hpd_file = filePath) # get total number of sampled trees
      
      hpdTrees = read.delim(file = filePath, comment.char = "#", skip = 3) # load hpd tree table
      
      outputDat[i, "nTrees"] = nrow(hpdTrees) # get number of unique trees in HPD interval 
      outputDat[i, "recovery"] = get_tree_recovery(hpdTrees$Tree, trueTree) # check if true tree is recovered in HPD interval 
      outputDat[i, "pRank1"] = hpdTrees[1, "Percent"] # check if true tree is recovered in HPD interval 
    }
    
    # print message
    print(paste0(tree, " file ", i, " done"))
  }
  
  outputDat = outputDat[order(outputDat$seed), ] # order according to simulation seed
  return(outputDat)
})

# store
names(hpdTreesOutput) = trees
saveRDS(hpdTreesOutput, file = paste0(outDir, "hpdTreesOutput.Rdat"))


# assess similarity between MCC and true tree topologies
mccOutput = lapply(trees, function(tree) {
  # find files with mcc tree
  infDir = paste0(dir, "inferenceOutput/")
  treesFiles = list.files(infDir, pattern = paste0(tree, "_[0-9]+.mcc.tree")) 
  n = length(treesFiles)
  
  # create output table
  outputDat = data.frame(matrix(nrow = n,  ncol = 4)) 
  colnames(outputDat) = c("seed", "wRF", "Nye", "shPI")
  
  for (i in 1:n) {
    # get mcc file and its simulation seed
    file = treesFiles[i]
    filePath = paste0(infDir, file)
    
    seed = as.numeric(str_extract(file, "[0-9]+"))
    outputDat[i, "seed"] = seed
    
    # get true tree
    trueTree = trueTrees[[tree]][[seed]]
    
    # get mcc tree
    mccTree = read.beast(file = filePath)
    mccTree = as.phylo(mccTree)
    
    # calculate tree similarity and distance metrics
    # weighted Robinson-Fould distance (considers branch lengths!)
    outputDat[i, "wRF"] = wRF.dist(mccTree, trueTree, normalize = T)
    # Generalized RFs: https://cran.r-project.org/web/packages/TreeDist/vignettes/Generalized-RF.html
    outputDat[i, "Nye"] = NyeSimilarity(mccTree, trueTree, normalize = T)  
    outputDat[i, "shPI"] = SharedPhylogeneticInfo(mccTree, trueTree, normalize = T)

    # print message
    print(paste0(tree, " file ", i, " done"))
  }
  
  outputDat = outputDat[order(outputDat$seed), ] # order according to simulation seed
  return(outputDat)
})

# store
names(mccOutput) = trees
saveRDS(mccOutput, file = paste0(outDir, "mccOutput.Rdat"))
