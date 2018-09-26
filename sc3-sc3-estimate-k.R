#!/usr/bin/env Rscript 
# Script to estimate the k size for SC3 clustering.

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
    help = "File name in which a processed SC3 object can be found."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store the SingleCellExperiment object with estimated k'."
  ),
  make_option(
    c("-t", "--output-text-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Text file name in which to store the estimated k'."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check if all parameters have been provided by the user

if ( ! file.exists(opt$input_object_file)){
  stop((paste('The SCE object', opt$input_object_file, 'does not exist')))
}

# Load common functions

suppressPackageStartupMessages(require(SC3))
suppressPackageStartupMessages(require(SingleCellExperiment))

# Read object

SingleCellExperiment <- readRDS(opt$input_object_file)

# Select K 

SingleCellExperiment <- sc3_estimate_k(SingleCellExperiment)

# Output to a serialized R object

saveRDS(SingleCellExperiment, file = opt$output_object_file)

# Output k to a text file

writeLines(as.character(metadata(SingleCellExperiment)$sc3$k_estimation), con = opt$output_text_file)
