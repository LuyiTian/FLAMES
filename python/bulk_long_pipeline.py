#!/usr/env python

import argparse
import os
import sys
import datetime
from parse_config import parse_json_config, print_config
from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions, remove_similar_tr, get_gene_flat, get_gene_blocks, group_bam2isoform
from gff3_to_fa import get_transcript_seq
from minimap2_align import minimap2_tr_align, gff3_to_bed12, minimap2_align, samtools_sort_index
from count_tr import parse_realigned_bam, parse_realigned_bam1, wrt_tr_to_csv, realigned_bam_coverage, parse_realigned_bam_raw
from filter_gff import annotate_filter_gff
from sc_long_pipeline import sc_long_pipeline
from merge_bulk_fq import merge_bulk_fq

__PROG = "FLAMES"
__AUTHOR = "Luyi Tian"
__VERSION = "0.1"
__MAN = \
    """
################################################################
# Program: {}
# Version {}
# Authors: {}
#
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
################################################################"""\
.format(__PROG, __VERSION, __AUTHOR)


def get_args():
    parser = argparse.ArgumentParser(
        prog=__PROG,
        description=__MAN,
        epilog="NOTE: make sure samtools is in PATH.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-a", "--gff3",
        help="The gene annotation in gff3 format.",
        type=str,
        required=True
    )

    parser.add_argument(
        "-i", "--fq_dir",
        help="folder containing fastq files.",
        type=str,
        default=True
    )

    parser.add_argument(
        "-b", "--inbam",
        help="aligned bam file (should be sorted and indexed). it will overwrite the `--infq` parameter and skip the first alignment step",
        type=str,
        default=""
    )

    parser.add_argument(
        "--outdir", "-o",
        help="directory to deposite all results in rootdir, use absolute path",
        type=str,
        required=True
    )

    parser.add_argument(
        "--genomefa", "-f",
        help="genome fasta file",
        type=str,
        required=True
    )

    parser.add_argument(
        "--minimap2_dir", "-m",
        help="directory contains minimap2, k8 and paftools.js program. k8 and paftools.js are used to convert gff3 to bed12.",
        type=str,
        required=False,
        default=""
    )

    parser.add_argument(
        "--config_file", "-c",
        help="json configuration files (default %(default)s)",
        type=str,
        default="config_sclr_nanopore_default.json"
    )

    parser.add_argument(
        "--downsample_ratio", "-d",
        help="downsampling ratio if performing downsampling analysis",
        type=float,
        default=1
    )

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = get_args()
    args.infq = "{}/merged.fastq.gz".format(args.outdir)

    print "Preprocessing bulk fastqs..."
    # creating output diretory for use by merge_bulk_fq
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        print("output directory not exist, create one:")
        print(args.outdir)
    merge_bulk_fq(args.fq_dir, "{}/barcode_anno.csv".format(args.outdir), args.infq)

    print "Running FLAMES pipeline..."
    sc_long_pipeline(args)
