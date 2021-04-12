# quantify transcript
import pysam
import os
from collections import Counter
import gzip
import numpy as np
import editdistance
from parse_gene_anno import parse_gff_tree


def umi_dedup(l, has_UMI):
    if has_UMI:
        l_cnt = Counter(l).most_common()
        if len(l_cnt)==1:
            return 1
        rm_umi = {}
        for ith in range(len(l_cnt)-1):
            for jth in range(len(l_cnt)-1,ith,-1):  # first assess the low abundant UMI
                if l_cnt[jth][0] not in rm_umi:
                    if editdistance.eval(l_cnt[ith][0],l_cnt[jth][0])<2:
                        rm_umi[l_cnt[jth][0]] = 1
        return(len(l_cnt)-len(rm_umi))
    else:
        return(len(l))


def wrt_tr_to_csv(bc_tr_count_dict, transcript_dict, csv_f, transcript_dict_ref=None, has_UMI=True):
    f = gzip.open(csv_f,"wb")
    all_tr = set()
    for bc in bc_tr_count_dict:
        all_tr.update(bc_tr_count_dict[bc].keys())
    all_tr = list(all_tr)
    f.write("transcript_id,gene_id,"+",".join([x for x in bc_tr_count_dict])+"\n" )
    tr_cnt = {}
    for tr in all_tr:
        cnt_l = [umi_dedup(bc_tr_count_dict[x][tr], has_UMI) if tr in bc_tr_count_dict[x] else 0 for x in bc_tr_count_dict ]
        tr_cnt[tr] = sum(cnt_l)
        if tr in transcript_dict:
            f.write("{},{},".format(tr,transcript_dict[tr].parent_id))
        elif (transcript_dict_ref is not None) and (tr in transcript_dict_ref):
            f.write("{},{},".format(tr,transcript_dict_ref[tr].parent_id))
        else:
            print "cannot find transcript in transcript_dict:", tr
            exit(1)
        f.write(",".join([str(x) for x in cnt_l])+"\n")
    f.close()
    return tr_cnt


def make_bc_dict(bc_anno):
    with open(bc_anno) as f:
        # skip header line
        f.readline()

        bc_dict = dict()
        for line in f:
            line_vals = line.rstrip().split(',')
            sample = line_vals[0]
            bc = line_vals[1]
            
            bc_dict[bc] = sample
        
    return(bc_dict)


def parse_realigned_bam(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage, **kwargs):
    """
    """
    fa_idx = dict((it.strip().split()[0],int(it.strip().split()[1]) ) for it in open(fa_idx_f))
    bc_tr_count_dict = {}
    bc_tr_badcov_count_dict = {}
    tr_cov_dict = {}
    read_dict = {}
    cnt_stat = Counter()
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    if "bc_file" in kwargs.keys() and kwargs["bc_file"] != "":
        bc_dict = make_bc_dict(kwargs["bc_file"])
    for rec in bamfile.fetch(until_eof=True):
        if rec.is_unmapped:
            cnt_stat["unmapped"] += 1
            continue
        map_st = rec.reference_start
        map_en = rec.reference_end
        tr = rec.reference_name
        tr_cov = float(map_en-map_st)/fa_idx[tr]
        tr_cov_dict.setdefault(tr,[]).append(tr_cov)
        if rec.query_name not in read_dict:
            read_dict.setdefault(rec.query_name,[]).append((tr, rec.get_tag("AS"), tr_cov, float(rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
        else:
            if rec.get_tag("AS")>read_dict[rec.query_name][0][1]:
                read_dict[rec.query_name].insert(0,(tr, rec.get_tag("AS"), tr_cov, float(rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
            elif rec.get_tag("AS")==read_dict[rec.query_name][0][1] and float(rec.query_alignment_length)/rec.infer_read_length()==read_dict[rec.query_name][0][3]:  # same aligned sequence
                if tr_cov>read_dict[rec.query_name][0][2]:  # choose the one with higher transcript coverage, might be internal TSS
                    read_dict[rec.query_name].insert(0,(tr, rec.get_tag("AS"), tr_cov, float(rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
            else:
                read_dict[rec.query_name].append((tr, rec.get_tag("AS"), tr_cov, float(rec.query_alignment_length)/rec.infer_read_length(), rec.mapping_quality))
        if tr not in fa_idx:
            cnt_stat["not_in_annotation"] += 1
            print tr, "not in annotation ???"
    tr_kept = dict((tr,tr) for tr in tr_cov_dict if len([it for it in tr_cov_dict[tr] if it > 0.9])>min_sup_reads)
    unique_tr_count = Counter(read_dict[r][0][0] for r in read_dict if read_dict[r][0][2]>0.9)
    for r in read_dict:
        tmp = read_dict[r]
        tmp = [it for it in tmp if it[0] in tr_kept]
        if len(tmp)>0:
            hit = tmp[0]  # transcript_id, pct_ref, pct_reads
        else:
            cnt_stat["no_good_match"] += 1
            continue
        bc, umi = r.split("#")[0].split("_")  # assume cleaned barcode
        if "bc_file" in kwargs.keys() and kwargs["bc_file"] != "":
            bc = bc_dict[bc]
        if len(tmp)==1 and tmp[0][4]>0:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
            cnt_stat["counted_reads"] += 1
        elif len(tmp)>1 and tmp[0][1]==tmp[1][1] and tmp[0][3]==tmp[1][3]:
            if hit[1] > 0.8:
                if bc not in bc_tr_count_dict:
                    bc_tr_count_dict[bc] = {}
                bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
                cnt_stat["counted_reads"] += 1
            else:
                cnt_stat["ambigious_reads"] += 1
                if bc not in bc_tr_badcov_count_dict:
                    bc_tr_badcov_count_dict[bc] = {}
                bc_tr_badcov_count_dict[bc].setdefault(hit[0], []).append(umi)
        elif hit[2] < min_tr_coverage or hit[3] < min_read_coverage:
            cnt_stat["not_enough_coverage"] += 1
            if bc not in bc_tr_badcov_count_dict:
                bc_tr_badcov_count_dict[bc] = {}
            bc_tr_badcov_count_dict[bc].setdefault(hit[0], []).append(umi)
        else:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(hit[0], []).append(umi)
            cnt_stat["counted_reads"] += 1
    print cnt_stat
    return bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept


def parse_realigned_bam1(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage):
    fa_idx = dict((it.strip().split()[0],int(it.strip().split()[1]) ) for it in open(fa_idx_f))
    bc_tr_count_dict = {}
    bc_tr_badcov_count_dict = {}
    tr_cov_dict = {}
    read_dict = {}
    cnt_stat = Counter()
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    for rec in bamfile.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_secondary: #or rec.mapping_quality==0:
            cnt_stat["not_counted"] += 1
            continue
        map_st = rec.reference_start
        map_en = rec.reference_end
        tr = rec.reference_name
        tr_cov = float(map_en-map_st)/fa_idx[tr]
        tr_cov_dict.setdefault(tr,[]).append(tr_cov)
        if tr not in fa_idx:
            cnt_stat["not_in_annotation"] += 1
            print tr, "not in annotation ???"
        bc, umi = rec.query_name.split("#")[0].split("_")  # assume cleaned barcode
        if float(map_en-map_st)/fa_idx[tr] < min_tr_coverage or float(rec.query_alignment_length)/rec.infer_read_length() < min_read_coverage:
            cnt_stat["not_enough_coverage"] += 1
            if bc not in bc_tr_badcov_count_dict:
                bc_tr_badcov_count_dict[bc] = {}
            bc_tr_badcov_count_dict[bc].setdefault(tr, []).append(umi)
        else:
            if bc not in bc_tr_count_dict:
                bc_tr_count_dict[bc] = {}
            bc_tr_count_dict[bc].setdefault(tr, []).append(umi)
        cnt_stat["counted_reads"] += 1
    tr_kept = dict((tr,tr) for tr in tr_cov_dict if len([it for it in tr_cov_dict[tr] if it > 0.9])>min_sup_reads)
    print cnt_stat
    return bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept


def parse_realigned_bam_raw(bam_in, fa_idx_f, min_sup_reads, min_tr_coverage, min_read_coverage):
    fa_idx = dict((it.strip().split()[0],int(it.strip().split()[1]) ) for it in open(fa_idx_f))
    bc_tr_count_dict = {}
    bc_tr_badcov_count_dict = {}
    tr_cov_dict = {}
    read_dict = {}
    cnt_stat = Counter()
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    for rec in bamfile.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_secondary: #or rec.mapping_quality==0:
            cnt_stat["not_counted"] += 1
            continue
        map_st = rec.reference_start
        map_en = rec.reference_end
        tr = rec.reference_name
        tr_cov = float(map_en-map_st)/fa_idx[tr]
        tr_cov_dict.setdefault(tr,[]).append(tr_cov)
        if tr not in fa_idx:
            cnt_stat["not_in_annotation"] += 1
            print tr, "not in annotation ???"
        bc, umi = rec.query_name.split("#")[0].split("_")  # assume cleaned barcode
        if bc not in bc_tr_count_dict:
            bc_tr_count_dict[bc] = {}
        bc_tr_count_dict[bc].setdefault(tr, []).append(umi)
        cnt_stat["counted_reads"] += 1
    tr_kept = dict((tr,tr) for tr in tr_cov_dict if len([it for it in tr_cov_dict[tr] if it > 0.9])>min_sup_reads)
    print cnt_stat
    return bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept


def tr_len_range(l):
    """
    [0,500] -> 0
    [500,1000] -> 1
    [1000,1500] -> 2
    [1500,2000] -> 3
    [2000,Inf] -> 4
    """
    return(min(4,l/500))


def realigned_bam_coverage(bam_in, fa_idx_f, coverage_dir):
    fa_idx = dict((it.strip().split()[0],int(it.strip().split()[1]) ) for it in open(fa_idx_f))
    left_clip_count = Counter()
    right_clip_count = Counter()
    tr_strand = Counter()
    bc_pct = {0:{},1:{},2:{},3:{},4:{}}
    bc_cov_pct = {0:[],1:[],2:[],3:[],4:[]}
    gene_pct = {0:[],1:[],2:[],3:[],4:[]}
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    for rec in bamfile.fetch(until_eof=True):
        if rec.is_unmapped or rec.is_supplementary or rec.is_secondary:
            continue
        bc, umi = rec.query_name.split("#")[0].split("_")  # assume cleaned barcode
        map_st = rec.reference_start
        map_en = rec.reference_end
        tr = rec.reference_name
        if float(map_en-map_st)/fa_idx[tr] < 0.3:
            continue
        if rec.cigar[0][0] == 4:  # BAM_CSOFT_CLIP
            left_clip_count[rec.cigar[0][1]] += 1
        if rec.cigar[-1][0] == 4:  # BAM_CSOFT_CLIP
            right_clip_count[rec.cigar[-1][1]] += 1
        tr_strand[rec.is_reverse] += 1
        if not rec.is_reverse:
            pass
        gene_pct[tr_len_range(fa_idx[tr])].append(float(map_en-map_st)/fa_idx[tr])
        bc_pct[tr_len_range(fa_idx[tr])].setdefault(bc, []).append(float(map_st-0)/fa_idx[tr])
        bc_pct[tr_len_range(fa_idx[tr])].setdefault(bc, []).append(float(map_en-0)/fa_idx[tr])
        bc_cov_pct[tr_len_range(fa_idx[tr])].append(float(map_en-map_st)/fa_idx[tr])
    print left_clip_count.most_common(30)
    print right_clip_count.most_common(30)
    print tr_strand
    print np.histogram(bc_pct[0][bc_pct[0].keys()[0]], bins=200, range=(0, 1))
    for i in bc_pct:
        coverage_f = open(os.path.join(coverage_dir,"transcript_cov_per_cell.{}.csv".format(i)), "w")
        for bc in bc_pct[i]:
            lhi, _ = np.histogram(bc_pct[i][bc], bins=200, range=(0, 1))
            coverage_f.write("{},".format(bc)+",".join(str(it) for it in lhi)+"\n")
        coverage_f.close()
    tr_cov_f = open(os.path.join(coverage_dir,"transcript_cov.csv"), "w")
    for i in gene_pct:
        lhi, _ = np.histogram(gene_pct[i], bins=200, range=(0, 1))
        tr_cov_f.write("{},".format(i)+",".join(str(it) for it in lhi)+"\n")
    tr_cov_f.close()



if __name__ == '__main__':
    #gff_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoforms/isoform_annotated.sample.nofilter.gff3"
    #csv_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoforms/transcript_count.sample.csv"
    #bam_in="/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/LT03_PromethION_10P.sample.realign.bam"
    """
    gff_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_test/isoform_annotated.gff3"
    csv_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_test/transcript_count.csv"
    bam_in = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_test/realign2transcript.bam"
    fa_idx_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_test/transcript_assembly.fa.fai"
    coverage_csv = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_test/per_cell_coverage.csv"
    tr_csv = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_test/per_tr_coverage.csv"
    
    gff_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all/isoform_annotated.gff3"
    csv_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all/transcript_count.csv"
    bam_in = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all/realign2transcript.bam"
    fa_idx_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all/transcript_assembly.fa.fai"
    coverage_csv = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all/per_cell_coverage.csv"
    tr_csv = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all/per_tr_coverage.csv"
    """
    data_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Rachel_scRNA_Aug19/isoform_out"
    bam_in = os.path.join(data_dir, "tmp_checkread.bam")
    fa_idx_f = os.path.join(data_dir,"transcript_assembly.fa.fai")
    bc_tr_count_dict, bc_tr_badcov_count_dict, tr_kept = parse_realigned_bam(bam_in, fa_idx_f, 3, 0.4, 0.4)
    print [len(bc_tr_count_dict[i]["ENSG00000277734.8_22493926_22552156_1"]) for i in bc_tr_count_dict if "ENSG00000277734.8_22493926_22552156_1" in bc_tr_count_dict[i]]
    print sum([len(bc_tr_count_dict[i]["ENSG00000277734.8_22493926_22552156_1"]) for i in bc_tr_count_dict if "ENSG00000277734.8_22493926_22552156_1" in bc_tr_count_dict[i]])



