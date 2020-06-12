import subprocess
import os
import pysam


def gff3_to_bed12(mm2_prog_path, gff3_file, bed12_file):
    if mm2_prog_path != "":
        cmd = "{_k8} {_paftools} gff2bed {_gff3} > {_bed}".format(
            _k8=os.path.join(mm2_prog_path, "k8"),
            _paftools=os.path.join(mm2_prog_path, "paftools.js"),
            _gff3=gff3_file,
            _bed=bed12_file)
    else:
        cmd = "paftools.js gff2bed {_gff3} > {_bed}".format(
            _gff3=gff3_file,
            _bed=bed12_file)
    print subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT)


def minimap2_align(mm2_prog_path, fa_file, fq_in, bam_out, no_flank=False, bed12_junc=None):
    """
    minimap2 align to genome
     -ax splice -t 12 -k14 $ref $read_fq | samtools view -bS -@ 4 -m 2G -o $out_bam -
samtools sort -@ 12 -o $sorted_bam $out_bam
samtools index $sorted_bam
    """
    if bed12_junc is not None:
        junc_cmd = "--junc-bed {} --junc-bonus 1".format(bed12_junc)
    else:
        junc_cmd = ""
    if no_flank:
        no_flank="--splice-flank=no"
    else:
        no_flank=""
    align_cmd = "{_prog} -ax splice -t 12 {_others} -k14 --secondary=no {_index} {_fq} | samtools view -bS -@ 4 -m 2G -o {_out} -  ".format(\
        _prog=os.path.join(mm2_prog_path, "minimap2"), _index=fa_file, _fq=fq_in, _out=bam_out, _others=" ".join([junc_cmd, no_flank]))
    print subprocess.check_output([align_cmd], shell=True, stderr=subprocess.STDOUT)


def samtools_sort_index(bam_in, bam_out):
    cmd = "samtools sort -@ 12 -o {_sorted_bam} {_in}".format(
        _sorted_bam=bam_out,
        _in=bam_in)
    print subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT)
    cmd = "samtools index {_sorted_bam}".format(_sorted_bam=bam_out)
    print subprocess.check_output([cmd], shell=True, stderr=subprocess.STDOUT)


def minimap2_tr_align(mm2_prog_path, fa_file, fq_in, bam_out):
    """
    minimap2 align to transcript
    """
    align_cmd = "{_prog} -ax map-ont -p 0.9 --end-bonus 10 -N 3 -t 12 {_index} {_fq} | samtools view -bS -@ 4 -m 2G -o {_out} -  ".format(\
        _prog=os.path.join(mm2_prog_path, "minimap2"), _index=fa_file, _fq=fq_in, _out=bam_out)
    # print align_cmd
    print subprocess.check_output([align_cmd], shell=True, stderr=subprocess.STDOUT)


if __name__ == '__main__':
    #mm2_prog_path = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/script/minimap2-2.11_x64-linux/"
    #genome_fa = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoforms/human_GRCh38_transcript.sample.fa"
    #fq_in="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/LT03_PromethION_10P.sorted.fastq.gz"
    #fq_in = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/LT03_PromethION_10P.sample.sort.fastq.gz"
    #bam_out = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/LT03_PromethION_10P.sample.realign.bam"
    #minimap2_tr_align(mm2_prog_path, genome_fa, fq_in, bam_out)
    mm2_prog_path="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/script/minimap2-2.17_x64-linux"
    fa_file="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCh38.primary_assembly.genome.fa"
    #fq_in="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/all_10p.fq"
    #tmp_bam="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/tmp_align_all.bam"
    #bam_out="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/align_all.sorted.bam"
    fq_in="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/MinION_revD/albacore_2.3.1/revd_fastq/sc_MinION_revD_pass_all.fq.gz"
    tmp_bam="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/MinION_revD/albacore_2.3.1/all_aligned.bam"
    bam_out="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/MinION_revD/albacore_2.3.1/all_aligned.sorted.bam"
    minimap2_align(mm2_prog_path, fa_file, fq_in, tmp_bam)
    samtools_sort_index(tmp_bam, bam_out)
    os.remove(tmp_bam)
