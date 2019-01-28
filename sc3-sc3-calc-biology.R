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
    c("-d", "--output-dir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "An output directory in which to write k-wise outputs."
  ),
  make_option(
    c("-k", "--ks"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "A comma-separated string or single value representing the number of clusters k to be used for SC3 clustering."
  ),
  make_option(
    c("-r", "--regime"),
    action = "store",
    default = NULL,
    type = 'character',
    help = "defines what biological analysis to perform. \"marker\" for marker genes, \"de\"\ for differentiall expressed genes and \"outl\" for outlier cells."
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name for R object of type 'SingleCellExperiment' from SC3 in which to store the consensus matrix."
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file', 'output_dir'))
saveRDS(opt, file = "opt.rds")

# Check parameter values defined
if ( ! file.exists(opt$input_object_file)){
  stop((paste('SC3 file object', opt$input_object_file, 'does not exist.')))
}


# Parse the ks to a vector
if ( is.null(opt$ks)){
  stop((paste('Please provide a k.')))
}else{
  ks <- wsc_parse_numeric(opt, 'ks')
}

# Once arguments are satisfcatory, load Scater package
suppressPackageStartupMessages(require(SC3))
suppressPackageStartupMessages(require(SingleCellExperiment))

# Input from serialized R object
SingleCellExperiment <- readRDS(opt$input_object_file)

# Calculate consensus matrix
SingleCellExperiment  <- sc3_calc_biology(object = SingleCellExperiment, ks = ks, regime = opt$regime)

# Get results by k
results_cols <- grep(opt$regime, colnames(rowData(SingleCellExperiment)), value = TRUE)

# Create output directory
dir.create(opt$output_dir)

for (k in ks){
  k_cols <- grep(paste0('_', k, '_'), results_cols, value = TRUE)
  results <- cbind(data.frame(cell = rownames(SingleCellExperiment)), rowData(SingleCellExperiment)[, k_cols, drop = FALSE])
  
  write.csv(results, file = file.path(opt$output_dir, paste(opt$regime, k, 'csv', sep='.')), row.names = FALSE, na='', quote = FALSE)
}

# Output to a serialized R object
saveRDS(SingleCellExperiment, file = opt$output_object_file)
