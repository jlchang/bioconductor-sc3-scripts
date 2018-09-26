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
    help = "File name in which a transformed SC3 'SingleCellExperiment' object has been stored after processed with sc3_calc_transfs()"
  ),
  make_option(
    c("-k", "--ks"),
    action = "store",
    default = NA,
    type = 'numeric',
    help = "A continuous range or a single integer representing the number of clusters - k - for SC3 clustering"
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for R object of type 'SingleCellExperiment' from SC3 in which to store the kmeans clustering as metadata."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values defined
if ( ! file.exists(opt$input_object_file)){
  stop((paste('SC3 file object', opt$input_object_file, 'does not exist.')))
}

if ( is.null(opt$k)){
  stop((paste('Please provide a k.')))
}

# Once arguments are satisfcatory, load Scater package
suppressPackageStartupMessages(require(SC3))

# Input from serialized R object
SingleCellExperiment <- readRDS(opt$input_object_file)

# calculate CPMs from raw count matrix
SingleCellExperiment  <- sc3_kmeans(object = SingleCellExperiment, k = opt$k)

# Output to a serialized R object
saveRDS(SingleCellExperiment, file = opt$output_object_file)
