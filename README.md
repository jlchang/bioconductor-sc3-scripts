# bioconductor-sc3-scripts [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/bioconductor-sc3-scripts/README.html)

In order to wrap SC3's internal workflow in any given workflow language, it's important to have scripts to call each of those steps, which is what this package provides

## Install

The recommended method for script installation is via a Bioconda recipe called bioconductor-sc3-scripts. 

With the [Bioconda channels](https://bioconda.github.io/#set-up-channels) configured the latest release version of the package can be installed via the regular conda install command:

```
conda install bioconductor-sc3-scripts
```

## Test installation

There is a test script included:

```
bioconductor-sc3-scripts-post-install-tests.sh
```

This downloads [a well-known test 10X dataset]('https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) and executes all of the scripts described below.

## Commands

The available wrapped SC3 functions are described below. Each script has usage instructions available via --help, consult function documentation in SC3 for further details. Wrappers are currently written for SC3 version 1.8.

### sc3_prepare(): prepare a SingleCellExperiment for SC3

See ?sc3_prepare for argument meanings.

```
sc3-sc3-prepare.R -i <SingleCellExperiment object as .rds> -f <gene_filter> -p <pct_dropout_min> \
    -q <pct_dropout_max> -d <d_region_min> -e <d_region_max> -n <svm_num_cells> -m <svm_max> \
    -t <n_cores> -s <rand_seed> -k <kmeans_nstart> -a <kmeans_iter_max> \
    -o <path to file where .rds output file will be stored>
```

### sc3_estimate_k(): Calculate k size for SC3 clustering

```
sc3-sc3-estimate-k.R -i <SC3 prepared SingleCellExperiment object as .rds> \
    -t <path to file where estimated k will be stored> \
    -o <path to file where .rds output file will be stored>
```

### sc3_calc_dists(): Calculate distances

```
sc3-sc3-calc-dists.R -i <SingleCellExperiment object from sc3_estimate_k() as .rds> \
    -o <path to file where .rds output file will be stored>
```

### sc3_calc_transfs():  Calculate transformations of the distance matrices

```
sc3-sc3-calc-transfs.R -i <SingleCellExperiment object from sc3_calc_dists() as .rds> \
    -o <path to file where .rds output file will be stored>
```

### sc3_kmeans(): Cluster transform matrix using k-means

```
sc3-sc3-kmeans.R -i <SingleCellExperiment object from sc3_calc_transfs() as .rds> \
    -k <k values to try, comma-separated> -o <path to file where .rds output file will be stored>
```

### sc3_calc_consens(): Calculate consensus clustering

```
sc3-calc-consens.R -i <SingleCellExperiment object from sc3_kmeans() as .rds> \
     -t <path to file where text output file will be stored> \
     -o <path to file where .rds output file will be stored>
```
