#!/usr/bin/env Rscript 

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))

# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

# parse options
option_list = list(
  make_option(
    c("-i", "--input-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which a SC3 'SingleCellExperiment' object has been stored after kmeans clustering."
  ),
  make_option(
    c("-t", "--output-text-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Text file name in which to store clusters, one column for every k value."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for R object of type 'SingleCellExperiment' from SC3 in which to store the consensus matrix."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_text_file'))

# Check parameter values defined
if ( ! file.exists(opt$input_object_file)){
  stop((paste('SC3 file object', opt$input_object_file, 'does not exist.')))
}

# Once arguments are satisfcatory, load Scater package
suppressPackageStartupMessages(require(SC3))
suppressPackageStartupMessages(require(SingleCellExperiment))

# Input from serialized R object
SingleCellExperiment <- readRDS(opt$input_object_file)

# Calculate consensus matrix
SingleCellExperiment  <- sc3_calc_consens(object = SingleCellExperiment)

# Derive the clusters matrix

clusters <- colData(SingleCellExperiment)[, grep('_clusters', colnames(colData(SingleCellExperiment))), drop = FALSE]

for (cluster_col in colnames(clusters)){
  k <- sub('.*_(\\d+)_.*', '\\1', cluster_col)
  
  # Write a summary for this K
  
  cluster_table <- data.frame(table(clusters[[cluster_col]]))
  colnames(cluster_table) <- c('Cluster', 'No. cells') 
  rownames(cluster_table) <- cluster_table$Cluster
  
  cat(
    paste(
      'Final cluster membership at k =', k,':\n'
    ),
    capture.output(cluster_table[, 2, drop = FALSE]),
    paste(length(which(is.na(clusters[[cluster_col]]))), 'cells fall outside clusters'),
    sep = '\n'
  )
}

# Output clusters to a tab-delimted file
clusters <- cbind(Cell=rownames(clusters), clusters)
write.table(clusters, file = opt$output_text_file, sep = "\t", row.names = FALSE, quote = FALSE, na='')

# Output to a serialized R object
saveRDS(SingleCellExperiment, file = opt$output_object_file)
