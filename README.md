# FLAMES
Full-length transcriptome splicing and mutation analysis

# Installation

The easiest way to install the dependencies for this package is using a conda environment.

```
conda create -n FLAMES python=2.7 samtools pysam minimap2 numpy editdistance -c bioconda -c conda-forge
```

The programs can then be run inside conda environment

```
git clone https://github.com/LuyiTian/FLAMES.git
conda activate FLAMES
FLAMES/python/sc_long_pipeline.py --help
FLAMES/python/bulk_long_pipeline.py --help
```
