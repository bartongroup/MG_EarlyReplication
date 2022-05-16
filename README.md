# The location and development of Replicon Cluster Domains in early replicating DNA

Software to accompany the manuscript de Costa et al.

## Usage

### On a Linux cluster

Copy all the FASTQ files from the archive to the `dna_seq/fastq` directory.

Create and activate a conda environment

```
cd dna_seq
conda env create -f env.yml
conda activate earlyrep_man
```

Run snakemake

```
./run_snake.sh
```

This will take a while - the fastq files are rather large.

### In RStudio

Once bedgraph files are created by snakemake, we suggest using RStudio. Start in the top project directory. The first step is to create environment using 'renv':

```
install.packages("renv")
renv::restore()
```

This will install all necessary packages. Run the `targets` pipeline.

```
targets::tar_make()
```

This will carry out all the calculations, create figures (as [targets](https://books.ropensci.org/targets/)) and output data in directory `tab`.

