#!/usr/bin/env python

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
        "-i", "--infq",
        help="input fastq file.",
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
        "--config_file","-c",
        help="json configuration files (default %(default)s)",
        type=str,
        default="config_sclr_nanopore_default.json"
        )
    parser.add_argument(
        "--downsample_ratio","-d",
        help="downsampling ratio if performing downsampling analysis",
        type=float,
        default=1
        )
    args = parser.parse_args()
    return args


def sc_long_pipeline(args):
    # parse configuration file

    if os.path.isfile(args.config_file):
        print("Use config file: {}".format(args.config_file))
        config_dict = parse_json_config(args.config_file)
    elif os.path.isfile(os.path.join(sys.path[0], args.config_file)):
        print("Use config file: {}".format(os.path.join(sys.path[0], args.config_file)))
        config_dict = parse_json_config(os.path.join(sys.path[0], args.config_file))
    else:
        print("Cannot find config file in current directory or script depository: {}".format(args.config_file))
        exit()
    print_config(config_dict)
    # check if files exist
    if args.downsample_ratio>1 or args.downsample_ratio<=0:
        print("downsample_ratio shoulw between 0 and 1: {}".format(args.downsample_ratio))
        exit()
    if not (os.path.isfile(args.infq) and os.path.isfile(args.gff3) and os.path.isfile(args.genomefa)):
        print("make sure all file exists:")
        print(args.infq)
        print(args.gff3)
        print(args.genomefa)
        exit()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
        print("output directory not exist, create one:")
        print(args.outdir)
    if args.inbam != "" and (not os.path.isfile(args.inbam)):
        print("make sure input inbam file exists:")
        print(args.inbam)
        exit()
    # output files:
    isoform_gff3 = os.path.join(args.outdir, "isoform_annotated.gff3")
    isoform_gff3_f = os.path.join(args.outdir, "isoform_annotated.filtered.gff3")
    FSM_anno_out = os.path.join(args.outdir, "isoform_FSM_annotation.csv")
    raw_splice_isoform = os.path.join(args.outdir, "splice_raw.gff3")
    tss_tes_stat = os.path.join(args.outdir, "tss_tes.bedgraph")
    transcript_fa = os.path.join(args.outdir, "transcript_assembly.fa")
    transcript_fa_idx = os.path.join(args.outdir, "transcript_assembly.fa.fai")
    tmp_bam = os.path.join(args.outdir, "tmp.align.bam")
    tmp_bed = os.path.join(args.outdir, "tmp.splice_anno.bed12")
    genome_bam = os.path.join(args.outdir, "align2genome.bam")
    realign_bam = os.path.join(args.outdir, "realign2transcript.bam")
    tr_cnt_csv = os.path.join(args.outdir, "transcript_count.csv.gz")
    tr_badcov_cnt_csv = os.path.join(args.outdir, "transcript_count.bad_coverage.csv.gz")
    print "Input parameters:"
    print "\tgene annotation:", args.gff3
    print "\tgenome fasta:", args.genomefa
    if args.inbam != "":
        print "\tinput bam:", args.inbam
        genome_bam = args.inbam
    else:
        print "\tinput fastq:", args.infq
    print "\toutput directory:", args.outdir
    print "\tdirectory contains minimap2:", args.minimap2_dir

    # align reads to genome
    if args.inbam == "" and config_dict["pipeline_parameters"]["do_genome_alignment"]:
        print "### align reads to genome using minimap2", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        if config_dict["alignment_parameters"]["use_junctions"]:
            gff3_to_bed12(args.minimap2_dir, args.gff3, tmp_bed)
        minimap2_align(args.minimap2_dir, args.genomefa, args.infq, tmp_bam, no_flank=config_dict["alignment_parameters"]["no_flank"], bed12_junc=tmp_bed if config_dict["alignment_parameters"]["use_junctions"] else None)
        samtools_sort_index(tmp_bam, genome_bam)
        os.remove(tmp_bam)
        if config_dict["alignment_parameters"]["use_junctions"]:
            os.remove(tmp_bed)
    else:
        print "### skip aligning reads to genome", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # find isoform
    print "### read gene annotation", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(args.gff3)
    transcript_to_junctions = {tr: blocks_to_junctions(transcript_to_exon[tr]) for tr in transcript_to_exon}
    remove_similar_tr(transcript_dict, gene_to_transcript, transcript_to_exon)
    gene_dict = get_gene_flat(gene_to_transcript, transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    if config_dict["pipeline_parameters"]["do_isoform_identification"]:
        print "### find isoforms", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        group_bam2isoform(genome_bam, isoform_gff3, tss_tes_stat, "", chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, args.genomefa,
        config=config_dict["isoform_parameters"], 
        downsample_ratio=args.downsample_ratio,
        raw_gff3=raw_splice_isoform if config_dict["global_parameters"]["generate_raw_isoform"] else None)
    else:
        print "### skip finding isoforms", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # get fasta
    #print "### generate transcript fasta file", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i = parse_gff_tree(isoform_gff3)
    ref_dict = {"chr_to_gene":chr_to_gene, "transcript_dict":transcript_dict, "gene_to_transcript":gene_to_transcript, "transcript_to_exon":transcript_to_exon}
    if not config_dict["realign_parameters"]["use_annotation"]:
        ref_dict = None
    get_transcript_seq(args.genomefa, transcript_fa, chr_to_gene_i, transcript_dict_i, gene_to_transcript_i, transcript_to_exon_i,ref_dict=ref_dict)

    # realign to transcript using minimap2
    if config_dict["pipeline_parameters"]["do_read_realignment"]:
        print "### realign to transcript using minimap2", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        minimap2_tr_align(args.minimap2_dir, transcript_fa, args.infq, tmp_bam)
        samtools_sort_index(tmp_bam, realign_bam)
        os.remove(tmp_bam)
    else:
        print "### skip read realignment", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # quantification
    if config_dict["pipeline_parameters"]["do_transcript_quantification"]:
        print "### generate transcript count matrix", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept = parse_realigned_bam(realign_bam, transcript_fa_idx, config_dict["isoform_parameters"]["Min_sup_cnt"], config_dict["transcript_counting"]["min_tr_coverage"], config_dict["transcript_counting"]["min_read_coverage"])
        #realigned_bam_coverage(realign_bam, transcript_fa_idx, args.outdir)
        tr_cnt = wrt_tr_to_csv(bc_tr_count_dict, transcript_dict_i, tr_cnt_csv, transcript_dict, config_dict["global_parameters"]["has_UMI"])
        wrt_tr_to_csv(bc_tr_badcov_count_dict, transcript_dict_i, tr_badcov_cnt_csv, transcript_dict, config_dict["global_parameters"]["has_UMI"])
        annotate_filter_gff(isoform_gff3,args.gff3,isoform_gff3_f,FSM_anno_out,tr_cnt,config_dict["isoform_parameters"]["Min_sup_cnt"])
    else:
        print "### skip transcript quantification", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
    args = get_args()
    sc_long_pipeline(args)
