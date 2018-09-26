#!/usr/bin/env Rscript 

# Run sc3_calc_dists() function. This function calculates distances between the
# cells. It creates and populates the following items of the sc3 slot of the
# metadata(object):
#
# distances - contains a list of distance matrices corresponding to Euclidean, Pearson and Spearman distances.

# Script to prepare SCE object for clustering with SC3.

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
    help = "File name in which a serialized R SingleCellExperiment object where object matrix found"
  ),
  make_option(
    c("-o", "--output-object-file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "File name in which to store serialized R object of type 'SingleCellExperiment'.'"
  )
)

opt <- wsc_parse_args(option_list, mandatory = c('input_object_file', 'output_object_file'))

# Check parameter values defined
if ( ! file.exists(opt$input_object_file)){
  stop((paste('File object or matrix', opt$input_object_file, 'does not exist')))
}

# Once arguments are satisfcatory, load SC3 package
suppressPackageStartupMessages(require(SC3))
suppressPackageStartupMessages(require(SingleCellExperiment))

# Input from serialized R object
SingleCellExperiment <- readRDS(opt$input_object_file)

# Run the function
# Seems to be necessary to convert sparse matrices to dense for this- see https://github.com/hemberg-lab/SC3/issues/53

assays(SingleCellExperiment) <- lapply(assays(SingleCellExperiment), function(x){
  if (! is.matrix(x)){
    as.matrix(x)
  }else{
    x
  }
})

SingleCellExperiment <- sc3_calc_dists(SingleCellExperiment)

# Output to a serialized R object

saveRDS(SingleCellExperiment, file = opt$output)
