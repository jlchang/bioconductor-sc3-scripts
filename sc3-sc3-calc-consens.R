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
    c("-d", "--output-cluster-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "An output directory in which to write clusters and marker genes."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for R object of type 'SingleCellExperiment' from SC3 in which to store the consensus matrix."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_cluster_dir'))

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

# Write a clusters file for each K

dir.create(opt$output_cluster_dir)

for (cluster_col in colnames(clusters)){
  k <- sub('.*(\\d+).*', '\\1', colnames(clusters))
  write.csv(data.frame(cell=rownames(clusters), cluster=clusters[[cluster_col]]), file = file.path(opt$output_cluster_dir, paste('clusters', k, 'csv', sep='.')), row.names = FALSE, na='', quote = FALSE)
}

# Output to a serialized R object
saveRDS(SingleCellExperiment, file = opt$output_object_file)
