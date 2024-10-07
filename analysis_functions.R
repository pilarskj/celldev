# ----------
# Functions for analysis_log.R script
# ----------

get_ntaxa = function(alignment_file) {
  "Get number of cells from simulated alignment (control)"
  
  # extract line from file
  line_4 = strsplit(readLines(alignment_file)[4], split = "=")
  
  # check that ntax is defined on that line
  if(grepl(line_4[[1]][1],pattern = "ntax")){
    ntaxa = as.numeric(strsplit(line_4[[1]][2], split = ";")[[1]][1])
    return(ntaxa)
  }else{
    stop("ntaxa not defined in line 4")
  }
}


load_alignment = function(alignment_file) {
  "Load simulated alignment from .nexus file 
  Create a data frame with rows - cells, columns - targets, values - edits"
  
  # read file 
  data = read.csv(alignment_file, skip = 11, header = F)
  
  # clean data frame
  data = data[-nrow(data), ] # remove last row
  data[nrow(data), ncol(data)] = str_replace(data[nrow(data), ncol(data)], ";", "") # remove ; in last cell
  rownames(data) = str_extract(sapply(str_split(data$V1, pattern = " "), "[", 1), "[0-9]+") # use cell labels as row names
  data$V1 = sapply(str_split(data$V1, pattern = " "), "[", 2) # remove cell labels from first column
  data$V1 = as.numeric(data$V1)
  data[, ncol(data)] = as.numeric(data[, ncol(data)])
  
  # sort according to cell labels
  data = data[order(as.numeric(rownames(data))), ] 
  return(data)
}


load_tapes = function(tape_files) {
  "Combine tapes from .nexus files to alignment"
  
  # load tapes
  data = sapply(tape_files, load_alignment, simplify = F)
  names(data) = str_extract(str_extract(names(data), "alignment_[0-9]+"), "[0-9]+")
  data = data[order(as.numeric(names(data)))]
  
  # rename target columns
  for (tape in names(data)) {
    colnames(data[[tape]]) = paste0("T", tape, "_", colnames(data[[tape]]))
  }
  
  # combine tapes
  data = bind_cols(data)
  return(data)
}


get_sumstat_names = function(parameters) {
  return(unlist(lapply(parameters, function(x){c(paste0(x,"_lower"), paste0(x, "_median"), paste0(x, "_upper"))})))
}


get_sumstats_from_log = function(log, parameters){
  "Extract inferred parameters from .log file and calculate summary statistics"
  
  # create output data frame
  output = data.frame(matrix(nrow = 1,  ncol = 3*length(parameters)))
  colnames(output) = get_sumstat_names(parameters)
  
  # for each parameter
  for (parameter in parameters){
    
    # get posterior estimates
    ix = match(parameter, colnames(log))
    posterior = log[, ix]
    
    # calculate median and HPD interval  
    hpd = hdi(posterior, credMass = 0.95)
    med = median(posterior)
    
    # save stats
    output[1, get_sumstat_names(parameter)] = c(hpd[[1]], med, hpd[[2]])
  }
  return(output)
}


get_ntrees = function(hpd_file) {
  "Get number of sampled trees in inference"
  
  # extract line from file
  line_total = strsplit(readLines(hpd_file)[3], split = ",")[[1]][2]
  
  # extract total number of trees from line
  ntrees = as.numeric(str_extract(line_total, "[0-9]+"))
  return(ntrees)
}


get_tree_recovery = function(hpdTrees, trueTree) {
  "Check if true tree is recovered in tree HPD interval"
  
  recovered = sapply(hpdTrees, function(x) {
    tree = read.tree(text = paste0(x, ";"))
    return(all.equal.phylo(target = trueTree, current = tree, use.edge.length = F))
  }, simplify = T, USE.NAMES = F)
  
  return(sum(recovered))
}


calc_hamming_dist <- function(x) {
  "Calculate pairwise Hamming distance between barcodes of an alignment"
  xt <- t(x)
  hdist = sapply(1:nrow(x), function(y) colSums(xt != xt[, y]))
  return(median(hdist, na.rm = T))
}


