# FLAMES

<img src="img/flames_logo.png" width="400">


Full-length transcriptome splicing and mutation analysis

# Installation

The easiest way to install the dependencies for this package is using a conda environment. The scripts can then be cloned from git.

```
conda create -n FLAMES \
    "python>=3.7" samtools pysam minimap2 numpy editdistance \
    -c bioconda -c conda-forge
git clone https://github.com/LuyiTian/FLAMES.git
```

# Usage for bulk data

Before using this software, remember to activate the FLAMES conda environment. The main scripts are the pipelines for single cells and bulk samples.

```
conda activate FLAMES
FLAMES/python/sc_long_pipeline.py --help
FLAMES/python/bulk_long_pipeline.py --help
```

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

# Usage for single cell data

For single cell, the first step is to find cell barcode in each long read. You can use a precompiled linux binary in `src/bin/match_cell_barcode`. If you encontered error, espically due to the different running environment. you can compile it from source code using the following command. Please use a C++ compiler that support C++11 features.

```
g++ -std=c++11 -lz -O2 -o match_cell_barcode ssw/ssw_cpp.cpp ssw/ssw.c match_cell_barcode.cpp kseq.h edit_dist.cpp
```

Then you can run `./match_cell_barcode` without argument to print the help message. `match_cell_barcode` requires a folder that contains all fastq file (can be `.gz` file), a file name/path for the statistics of barcode matching, a csv file of cell barcode that used as reference and a file name/path to the output fastq.gz file. The cell barcode for 10x will be in `filtered_feature_bc_matrix/barcodes.tsv.gz`, please unzip the file and use `barcodes.tsv` as input. It also support scPipe cell barcode annotation file generated from `sc_detect_bc` function. In `match_cell_barcode`, the flanking sequence (`CTACACGACGCTCTTCCGATCT`) is aligned to the first 30000 reads to identify the regions where cell barcode is likely to be found within. Next, sequences within this region are matched to barcodes in `barcodes.tsv`, allowing `MAX_DIST` (the 5th argument of `match_cell_barcode`) hamming distances. Reads that are successfully matched with a barcode are reported in the `barcode hm match` count. Reads that could not be matched in the previous step are aligned to the flanking sequence to identify the location of barcode individually, and barcode matching is done with up to `MAX_DIST` levenshtein distances (allowing indels). Reads that are matched by this step is reported by the `fuzzy match` count.

Next, after you get the fastq file from `match_cell_barcode`. you could run `sc_long_pipeline.py` with the following command.


```
usage: FLTSA [-h] -a GFF3 [-i INFQ] [-b INBAM] --outdir OUTDIR --genomefa
             GENOMEFA --minimap2_dir MINIMAP2_DIR [--config_file CONFIG_FILE]
             [--downsample_ratio DOWNSAMPLE_RATIO]

# semi-supervised isoform detection and annotation from long read data.
# output:
# outdir:
#   transcript_count.csv.gz   // transcript count matrix
#   isoform_annotated.filtered.gff3 // isoforms in gff3 format
#   transcript_assembly.fa // transcript sequence from the isoforms
#   align2genome.bam       // sorted bam file with reads aligned to genome
#   realign2transcript.bam // sorted realigned bam file using the
#                            transcript_assembly.fa as reference
#   tss_tes.bedgraph       // TSS TES enrichment for all reads (for QC)
################################################################

optional arguments:
  -h, --help            show this help message and exit
  -a GFF3, --gff3 GFF3  The gene annotation in gff3 format.
  -i INFQ, --infq INFQ  input fastq file.
  -b INBAM, --inbam INBAM
                        aligned bam file (should be sorted and indexed). it
                        will overwrite the `--infq` parameter and skip the
                        first alignment step
  --outdir OUTDIR, -o OUTDIR
                        directory to deposite all results in rootdir, use
                        absolute path
  --genomefa GENOMEFA, -f GENOMEFA
                        genome fasta file
  --minimap2_dir MINIMAP2_DIR, -m MINIMAP2_DIR
                        directory contains minimap2, k8 and paftools.js
                        program. k8 and paftools.js are used to convert gff3
                        to bed12.
  --config_file CONFIG_FILE, -c CONFIG_FILE
                        json configuration files (default
                        config_sclr_nanopore_default.json)
  --downsample_ratio DOWNSAMPLE_RATIO, -d DOWNSAMPLE_RATIO
                        downsampling ratio if performing downsampling analysis
```

# configuration file

FLAMES provides a default set of parameters, but can be changed by the configuration JSON file.
The `pipeline_parameters` section specifies which step to be excuated in the pipeline, by default you should go though all steps.
The `isoform_parameters` section determines the results isoform detection, some key parameters include `Min_sup_cnt` which means transcript with less read aligned than `Min_sup_cnt` will be discarded, `MAX_TS_DIST` which will merge transcripts with the same intron chain and TSS/TES distance less than `MAX_TS_DIST`. `strand_specific` will specify whether the the read is strand specific, such as the reads are in the same strand as the mRNA (1) or the reverse complement (-1), or the reads are not strand specific (0), which means the method will determine the strand information based on reference annotation.
