# ----------
# Script for the analysis of simulated alignments from .nexus files and
# BEAST inference results from .log, .hpd.trees and .mcc.tree files

# During inference, the following parameters were estimated from the simulated data: 
# birth rate, death rate, migration rates, edit rate, tree length, tree height
# The sampling proportion, the origin, the time interval of scarring and the scarring probabilities were fixed.
# ----------

# load packages
library(tracerer)
library(HDInterval)
library(ape)
library(treeio)
library(TreeDist) 
library(phangorn)
library(TreeTools)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

rm(list = ls()) # clear environment

# adapt paths
code_dir = '~/Projects/celldev'
data_dir = '~/Projects/celldev_data'
setwd(code_dir)

# load helper functions
source("analysis_functions.R")


### SCENARIO
nSim = 25

# lineage recorder
method = "Typewriter"

# multitype model prototype
prototype = "distinct"

# inferred parameters
parameters = c("treeHeight", "treeLength", "editRate", 
               "birthRate.t0", "birthRate.t1", "birthRate.t2", "birthRate.t3",
               "deathRate.t0", "deathRate.t1", "deathRate.t2", "deathRate.t3", 
               "migrationRate.t1", "migrationRate.t2", "migrationRate.t3")
nParams = length(parameters)

# tree parameters
treeParams = read.csv(paste0("heterog/simulation_trees/tree_params_", prototype, ".csv"))  

# directory for current setup
dir = paste0(data_dir, "/heterog/", method, "/", prototype)

# output directory
outDir = file.path(dir, "analysisOutput")
if (!dir.exists(outDir)) {dir.create(outDir)}



### SIMULATED ALIGNMENTS
simDir = file.path(dir, "simulationOutput")
simFiles = list.files(simDir, pattern = "*.nexus", full.names = T)
n = length(simFiles)

# load alignments
if (method == "TiDe"){
  alignments = sapply(simFiles, load_alignment, simplify = F)
  names(alignments) = str_extract(names(alignments), "[0-9]+") # extract simulation seed
} else if (method == "Typewriter"){
  seeds = unique(str_extract(sapply(str_split(simFiles, "[.]"), "[[", 1), "[0-9]+"))
  alignments = sapply(seeds, function(i) {
    tape_files = simFiles[str_detect(simFiles, paste0("_", i, "[.]alignment"))]
    data = load_tapes(tape_files)
    return(data)
  }, simplify = F, USE.NAMES = T)
}

# get unique barcodes per alignment
barcodes = lapply(alignments, function(alignment) {
  alignment %>% group_by_all() %>% summarize(count=n()) %>% arrange(desc(count)) %>%
    suppressMessages() })

# calculate statistics per alignment
simStats = data.frame(matrix(nrow = min(n, nSim), ncol = 6))
colnames(simStats) = c("seed", "nEdited", "nEdits", "nBarcodes", "nCells", "barcodeDiv")

simStats$seed = names(alignments)

# get median number of edited targets/tapes per cell
if (method == "TiDe"){
  simStats$nEdited = sapply(alignments, function(a) {median(rowSums(a != 0)) })
} else if (method == "Typewriter"){
  simStats$nEdited = sapply(alignments, function(a) {median(rowSums(a %>% select(contains("V1")) != 0)) })
}

# get number of different edits made across targets and cells
simStats$nEdits = sapply(alignments, function(a) {length(unique(unlist(a)[unlist(a) != 0])) })

# get number of unique barcodes
simStats$nBarcodes = sapply(barcodes, nrow)

# check if number of cells in simulations corresponds to number of tips in tree
if (method == "TiDe"){
  simStats$nCells = sapply(simFiles, get_ntaxa)
} else if (method == "Typewriter"){
  simStats$nCells = sapply(seeds, function(i) {
    get_ntaxa(simFiles[str_detect(simFiles, paste0("_", i, "[.]alignment_1[.]"))])})
}
print(all(simStats$nCells == treeParams[match(simStats$seed, treeParams$seed), "ntips"]))

# calculate barcode diversity: number of barcodes per cell (how many cells share a barcode?)
simStats$barcodeDiv = simStats$nBarcodes / simStats$nCells

# sort by seed
simStats = simStats %>% arrange(as.numeric(seed))

saveRDS(simStats, file = file.path(outDir, "simStats.Rdat"))



### INFERENCE RESULTS
# summarize inference output for each tree
infDir = file.path(dir, "inferenceOutput")
logFiles = list.files(infDir, pattern = "*.inference.log") # find log files
n = length(logFiles)
  
infOutput = data.frame(matrix(nrow = n,  ncol = 3*nParams + 3)) # create output table
colnames(infOutput) = c(get_sumstat_names(parameters), "seed", "minESS", "chainLength")
  
for (i in 1:n){
  # get log file and its simulation seed
  file = logFiles[i]
  infOutput[i, "seed"] = str_extract(file, "[0-9]+")
    
  # parse beast log 
  log = remove_burn_ins(parse_beast_tracelog_file(file.path(infDir, file)), burn_in_fraction = 0.1)
  
  # rename relevant parameters
  log = log %>% rename(treeHeight = tree.height, treeLength = tree.treeLength,
                       birthRate.t0 = birthRate.1, birthRate.t1 = birthRate.2, birthRate.t2 = birthRate.3, birthRate.t3 = birthRate.4,
                       deathRate.t0 = deathRate.1, deathRate.t1 = deathRate.2, deathRate.t2 = deathRate.3, deathRate.t3 = deathRate.4) %>%
    select(where(~ any(. != 0)))
  
  if (prototype == "distinct") {
    log = log %>% rename(migrationRate.t1 = migrationRate.distinct.1,
                         migrationRate.t2 = migrationRate.distinct.2,
                         migrationRate.t3 = migrationRate.distinct.3)
  } else if (prototype == "hierarchical") {
    log = log %>% rename(migrationRate.t1 = migrationRate.hierarchical.1,
                         migrationRate.t2 = migrationRate.hierarchical.5,
                         migrationRate.t3 = migrationRate.hierarchical.9)
  }
  
  # check ESS
  ess = calc_esses(log, sample_interval = 5000)
  infOutput[i, "minESS"] = min(ess)
  
  # get chain length
  infOutput[i, "chainLength"] = log$Sample[nrow(log)]
  
  # save summary stats in output data frame
  infOutput[i, 1:(3*nParams)] = get_sumstats_from_log(log, parameters)
    
  # print message
  print(paste0("file ", i, " done"))
}
  
infOutput$seed = as.numeric(infOutput$seed)
infOutput = infOutput[order(infOutput$seed), ] # order according to simulation seed

saveRDS(infOutput, file = file.path(outDir, "infOutput.Rdat"))

# infOutput = readRDS(file.path(outDir, "infOutput.Rdat"))


### INFERENCE PERFORMANCE
infOutput = infOutput %>% filter(minESS >= 200) # filter valid simulations
print(nrow(infOutput))

# true parameters
true_parameters = c("editRate" = 0.05, 
                    "birthRate.t0" = 0.2, "birthRate.t1" = 0.18, "birthRate.t2" = 0.18, "birthRate.t3" = 0.18,
                    "deathRate.t0" = 0.01, "deathRate.t1" = 0.03, "deathRate.t2" = 0.03, "deathRate.t3" = 0.03)
if (prototype == "distinct") {
  true_parameters = c(true_parameters, "migrationRate.t1" = 0.04, "migrationRate.t2" = 0.03, "migrationRate.t3" = 0.02)
} else if (prototype == "hierarchical") {
  true_parameters = c(true_parameters, "migrationRate.t1" = 0.08, "migrationRate.t2" = 0.06, "migrationRate.t3" = 0.04)
}

# calculate width of priors
# priors = list()
# priors$editRate = c(hdi(qlnorm, 0.95, meanlog = -3, sdlog = 1),
#                     "median" = qlnorm(0.5, meanlog = -3, sdlog = 1))
# priors$birthRate = c(hdi(qlnorm, 0.95, meanlog = -1.5, sdlog = 1),
#                      "median" = qlnorm(0.5, meanlog = -1.5, sdlog = 1))
# priors$deathRate = c(hdi(qexp, 0.95, rate = 25),
#                      "median" = qexp(0.5, rate = 25))
# priors$migrationRate = c(hdi(qlnorm, 0.95, meanlog = -3.25, sdlog = 1),
#                          "median" = qlnorm(0.5, meanlog = -3.25, sdlog = 1))
# priors_hpd_width = sapply(priors, function(x) {as.numeric(x["upper"] - x["lower"])})
priors_hpd_width = c("editRate" = 0.2569439, 
                     "birthRate.t0" = 1.1515428, "birthRate.t1" = 1.1515428, "birthRate.t2" = 1.1515428, "birthRate.t3" = 1.1515428,
                     "deathRate.t0" = 0.1198293, "deathRate.t1" = 0.1198293, "deathRate.t2" = 0.1198293, "deathRate.t3" = 0.1198293,
                     "migrationRate.t1" = 0.2001081, "migrationRate.t2" = 0.2001081, "migrationRate.t3" = 0.2001081)

# specify performance metrics
perf_metrics = c("coverage", "abs_error_rel", "bias_rel", "hpd_width_rel", "hpd_proportion")

infPerformance = sapply(perf_metrics, function(x) {
  perf = data.frame(matrix(nrow = nrow(infOutput), ncol = length(parameters) + 1, data = NA))
  colnames(perf) = c("seed", parameters)
  perf$seed = infOutput$seed 
  return(perf)}, simplify = F, USE.NAMES = T)
  
for (i in 1:nrow(infOutput)) {
  for(parameter in parameters) {
    seed = infOutput[i, "seed"]
    
    # get true parameter
    if (parameter %in% c("treeHeight", "treeLength")) {
      truth = treeParams[treeParams$seed == seed, parameter]
    } else {
      truth = true_parameters[parameter]
    }
    
    # get estimated parameter
    median = infOutput[i, paste0(parameter, "_median")]
    hpd_lower = infOutput[i, paste0(parameter, "_lower")]
    hpd_upper = infOutput[i, paste0(parameter, "_upper")]
    
    # calculate normalized performance metrics
    infPerformance$coverage[i, parameter] = ifelse(truth >= hpd_lower & truth <= hpd_upper, TRUE, FALSE)
    infPerformance$abs_error_rel[i, parameter] = abs(median - truth) / truth
    infPerformance$bias_rel[i, parameter] = (median - truth) / truth
    infPerformance$hpd_width_rel[i, parameter] = (hpd_upper - hpd_lower) / truth
    
    # calculate information gain
    infPerformance$hpd_proportion[i, parameter] = (hpd_upper - hpd_lower) / priors_hpd_width[parameter]
  }
}
    
saveRDS(infPerformance, file = file.path(outDir, "infPerformance.Rdat"))


### TREE INFERENCE
# load true tree topologies from .newick files
trueTreesFiles = list.files(file.path(data_dir, paste0("Trees/", prototype)), 
                            pattern = ".typed.node.trees", full.names = T)
trueTrees = sapply(trueTreesFiles, function(tree) {as.phylo(read.beast(tree))}, simplify = F)
names(trueTrees) = str_extract(names(trueTrees), "[0-9]+")
trueTypes = sapply(trueTreesFiles, function(tree) {get.data(read.beast(tree)) %>% arrange(node)}, simplify = F)
names(trueTypes) = str_extract(names(trueTypes), "[0-9]+")


# analyse 95% credible sets of trees
hpdTreesFiles = list.files(infDir, pattern = ".hpd.trees") # find files with tree HPDs
n = length(hpdTreesFiles)
   
hpdTreesOutput = data.frame(matrix(nrow = n,  ncol = 4)) # create output table
colnames(hpdTreesOutput) = c("seed", "nTrees", "nTotal", "recovery")

for (i in 1:n) {
  # get trees file and its simulation seed
  file = hpdTreesFiles[i]
  filePath = file.path(infDir, file)
  
  seed = str_extract(file, "[0-9]+")
  hpdTreesOutput[i, "seed"] = seed
  
  trueTree = trueTrees[[seed]]
  
  # if file not empty
  if (!is.na(readLines(filePath)[4])) {
    hpdTreesOutput[i, "nTotal"] = get_ntrees(hpd_file = filePath) # get total number of sampled trees
    
    hpdTrees = read.delim(file = filePath, comment.char = "#", skip = 3) # load hpd tree table
    
    hpdTreesOutput[i, "nTrees"] = nrow(hpdTrees) # get number of unique trees in HPD interval
    hpdTreesOutput[i, "recovery"] = get_tree_recovery(hpdTrees$Tree, trueTree) # check if true tree is recovered in HPD interval
  }
  
  # print message
  print(paste0("file ", i, " done"))
}

hpdTreesOutput$seed = as.numeric(hpdTreesOutput$seed)
hpdTreesOutput = hpdTreesOutput[order(hpdTreesOutput$seed), ] # order according to simulation seed
saveRDS(hpdTreesOutput, file = file.path(outDir, "hpdTreesOutput.Rdat"))


# assess similarity between MCC and true tree topologies
mccTreesFiles = list.files(infDir, pattern = "*.mcc.tree") # find files with mcc tree
n = length(mccTreesFiles)

mccOutput = data.frame(matrix(nrow = n,  ncol = 5)) # create output table
colnames(mccOutput) = c("seed", "wRF", "Nye", "shPI", "prop_type_correct")

for (i in 1:n) {
  # get mcc file and its simulation seed
  file = mccTreesFiles[i]
  filePath = file.path(infDir, file)
  
  seed = str_extract(file, "[0-9]+")
  mccOutput[i, "seed"] = seed
  
  # get true tree and types
  trueTree = trueTrees[[seed]]
  trueType = trueTypes[[seed]]
  
  # get mcc tree and types
  mccTree = read.beast(file = filePath)
  mccType = get.data(mccTree) %>% select(type, node) %>% arrange(node)
  mccTree = as.phylo(mccTree)
  
  # calculate tree similarity and distance metrics
  # weighted Robinson-Fould distance (considers branch lengths!)
  mccOutput[i, "wRF"] = wRF.dist(mccTree, trueTree, normalize = T)
  # generalized RF metrics
  mccOutput[i, "Nye"] = NyeSimilarity(mccTree, trueTree, normalize = T)  
  mccOutput[i, "shPI"] = SharedPhylogeneticInfo(mccTree, trueTree, normalize = T)
  
  # calculate proportion of types correctly assigned to internal nodes in mcc tree
  compareTrees = comparePhylo(mccTree, trueTree)
  commonNodes = compareTrees$BT %>%
    mutate(mccTree = str_extract(str_extract(mccTree, "\\([0-9]+\\)"), "[0-9]+"),
           trueTree = str_extract(str_extract(trueTree, "\\([0-9]+\\)"), "[0-9]+"))
  
  mccOutput[i, "prop_type_correct"] = sum(mccType[commonNodes$mccTree, "type"] == trueType[commonNodes$trueTree, "type"]) / nrow(commonNodes)
  
  # print message
  print(paste0("file ", i, " done"))
}

mccOutput$seed = as.numeric(mccOutput$seed)
mccOutput = mccOutput[order(mccOutput$seed), ] # order according to simulation seed
saveRDS(mccOutput, file = file.path(outDir, "mccOutput.Rdat"))
