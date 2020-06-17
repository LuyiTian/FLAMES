# FLAMES

<img src="img/flames_logo.png" width="400">


Full-length transcriptome splicing and mutation analysis

# Installation

The easiest way to install the dependencies for this package is using a conda environment. The scripts can then be cloned from git.

```
conda create -n FLAMES \
    python=2.7 samtools pysam minimap2 numpy editdistance \
    -c bioconda -c conda-forge
git clone https://github.com/LuyiTian/FLAMES.git
```

# Usage

Before using this software, remember to activate the FLAMES conda environment. The main scripts are the pipelines for single cells and bulk samples.

```
conda activate FLAMES
FLAMES/python/sc_long_pipeline.py --help
FLAMES/python/bulk_long_pipeline.py --help
```

# Example

An example has been included with a small subset of SIRV data in the examples folder.

```
PATH=$PWD/FLAMES/python:$PATH

cd examples/SIRV
bulk_long_pipeline.py \
    --gff3 data/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf \
    --genomefa data/SIRV_isoforms_multi-fasta_170612a.fasta \
    --outdir FLAMES_output \
    --config_file data/SIRV_config.json \
    --fq_dir data/fastq

# output data is in FLAMES_output
ls FLAMES_output
```
