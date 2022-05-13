# The location and development of Replicon Cluster Domains in early replicating DNA

Software to accompany the manuscript.

## Usage

### On a Linux cluster

First, create and activate a conda environment:

```
cd dna_seq
conda env create -f env.yml
conda activate earlyrep_man
```

To run snakemake

```
./run_snake.sh
```

This will take a while - the fastq files are rather large.

### In RStudio

Once bedgraph files are created by snakemake, we suggest using RStudio. Start it in the top project directory. The first step is to create environment using 'renv':

```
install.packages("renv")
renv::restore()
```

Once all the packages are installed, run the `targets` pipeline.

```
targets::tar_make()
```

This will carry out all the calculations, create figures and output data.

