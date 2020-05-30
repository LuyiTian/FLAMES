# detect mutations in bam files.
import os
import pysam
import gzip
import numpy as np
from scipy.stats import hypergeom
from collections import Counter, namedtuple
from parse_gene_anno import parse_gff_tree
from sc_longread import blocks_to_junctions, get_gene_flat, get_gene_blocks


def get_fa(fn):
    ch = ""
    seq = []
    for line in open(fn):
        if line[0] == ">":
            if ch != "":
                yield ch, "".join(seq)
            ch = line[1:].strip().split()[0]
            seq = []
        else:
            seq.append(line.strip().upper())
    yield ch, "".join(seq)


def seq_entropy(seq):
    res = 0.
    for st in list(set(seq)):
        p = float(seq.count(st))/len(seq)
        res += -p*np.log(p)
    return res


def find_homo_regions(fa_seq, chr_bl, min_len=3, min_gap=1):
    """
    find regions with at least `min_len` homopolymers and 
    call +/- `min_gap` in surrounding regions as homo-regions as well.
    """
    homo_dict = {}
    for bl in chr_bl:
        i=bl.s
        while i < bl.e-min_len-1:
            if fa_seq[i:(i+min_len)] == "A" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="A":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "A"
                i = j
            elif fa_seq[i:(i+min_len)] == "T" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="T":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "T"
                i = j
            elif fa_seq[i:(i+min_len)] == "G" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="G":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "G"
                i = j
            elif fa_seq[i:(i+min_len)] == "C" * min_len:
                j= i+min_len
                while j<bl.e-1 and fa_seq[j]=="C":
                    j += 1
                for ix in range(max(0,i-min_gap),min(j+min_gap,bl.e-1)):
                    homo_dict[ix] = "C"
                i = j
            else:
                i += 1
    return homo_dict


def update_corr_cnt(int_l, cb_corr_cnt):
    for i in range(len(int_l)-1):
        for j in range(i+1, len(int_l)):
            cb_corr_cnt[(int_l[i],int_l[j])] += 1
            cb_corr_cnt[(int_l[j],int_l[i])] += 1




def realigned_bam_allele_coverage(bam_in, chr_to_blocks, fa_f, cov_bin_f, cb_seq_dict, vcf_f=None, min_cnt=150,min_cov=100,report_pct=(0.1,0.9) ):
    c2i = {"A":0, "C":1, "G":2, "T":3}  # four array.arrays of the same length in order A C G T
    fa_dict={}
    vcf_dict={}
    acc_pct = []
    cb_seq_set = set(cb_seq_dict.keys())
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    #vcf_in = pysam.VariantFile(vcf_f)
    cb_corr_cnt = Counter()
    vcf_c = 0
    vcf_not_c = 0
    for ch in chr_to_blocks:
        print ch
        #if ch != "chr15":
        #    continue
        homo_dict = find_homo_regions(fa_dict[ch], chr_to_blocks[ch])
        for ith, bl in enumerate(chr_to_blocks[ch]):
            cnt = bamfile.count(ch, bl.s, bl.e)
            #try:
            #    vcf_dict = dict((it.pos-1, it) for it in vcf_in.fetch(ch[3:], bl.s, bl.e))
            #except:
            #    print ch[3:], "not in vcf.  ",ch
            if cnt < min_cnt:
                continue
            acc_pct_tr = []
            cov = bamfile.count_coverage(ch, bl.s, bl.e,
            quality_threshold=0)  # four array.arrays of the same length in order A C G T
            for i in range(10, len(cov[0])-10):  # ignore the bases at the beginning and the end
                tot =  float(cov[0][i]+cov[1][i]+cov[2][i]+cov[3][i])
                if (bl.s+i not in homo_dict) and tot>min_cov:
                    acc_pct_tr.append((bl.s+i, cov[c2i[fa_dict[ch][bl.s+i]]][i]/tot, [("A",cov[0][i]),("C",cov[1][i]),("G",cov[2][i]),("T",cov[3][i])] ))
            if len(acc_pct_tr)>10:
                for ix, pct in enumerate(acc_pct_tr):
                    if ix > 1 and ix < len(acc_pct_tr)-1:
                        if report_pct[0]<pct[1]<report_pct[1]:
                            if acc_pct_tr[ix-1][1]>0.95 and acc_pct_tr[ix+1][1]>0.95:
                                if seq_entropy(fa_dict[ch][(pct[0]-10):pct[0]])<1 or seq_entropy(fa_dict[ch][pct[0]:(pct[0]+10)])<1:  # ignore homo regions
                                    continue
                                tmp_atcg_set = {}
                                tmp_set = Counter()
                                for pileupcolumn in bamfile.pileup(ch, pct[0], pct[0]+1,truncate=True, min_base_quality=0,ignore_overlaps=False):
                                    c_keep = 0
                                    c_del = 0
                                    for pileupread in pileupcolumn.pileups:
                                        if not pileupread.is_del and not pileupread.is_refskip:
                                            c_keep += 1
                                            cb_seq, umi_seq = pileupread.alignment.query_name.split("#")[0].split("_")
                                            if cb_seq in cb_seq_set:
                                                tmp_atcg_set.setdefault(pileupread.alignment.query_sequence[pileupread.query_position],Counter())[cb_seq] += 1
                                                tmp_set[cb_seq] += 1
                                        else:
                                            c_del += 1
                                if c_keep/float(c_keep+c_del)<0.7:
                                    continue
                                bs = tmp_atcg_set.keys()
                                for b in bs:
                                    tmp_atcg_set[b] = set(it for it in tmp_atcg_set[b] if tmp_atcg_set[b][it]>1)
                                tmp_set = set(it for it in tmp_set if tmp_set[it]>1)
                                pct[2].sort(key=lambda x:x[1],reverse=True)
                                lead_b = pct[2][0][0]  # only look at most enriched two possibilities
                                snd_b = pct[2][1][0]
                                if not (len(tmp_atcg_set[lead_b])>10 and len(tmp_atcg_set[snd_b])>10):
                                    continue
                                rv = hypergeom(len(tmp_set), len(tmp_atcg_set[lead_b]), len(tmp_atcg_set[snd_b]))
                                if rv.pmf(len(tmp_atcg_set[lead_b].intersection(tmp_atcg_set[snd_b])))<0.000001 and len(tmp_atcg_set[lead_b].intersection(tmp_atcg_set[snd_b]))<0.9*min(len(tmp_atcg_set[lead_b]), len(tmp_atcg_set[snd_b])):
                                    #print ("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
                                    tmp_set = tmp_atcg_set[lead_b] - tmp_atcg_set[snd_b]  # x not y
                                    if len(tmp_set)>1:
                                        update_corr_cnt(list(tmp_set), cb_corr_cnt)
                                    tmp_set = tmp_atcg_set[snd_b] - tmp_atcg_set[lead_b]  # y not x
                                    if len(tmp_set)>1:
                                        update_corr_cnt(list(tmp_set), cb_corr_cnt)
                                    tmp_set = tmp_atcg_set[lead_b] & tmp_atcg_set[snd_b]  # x and y
                                    if len(tmp_set)>1:
                                        update_corr_cnt(list(tmp_set), cb_corr_cnt)
                                    if pct[0] in vcf_dict:
                                        vcf_c += 1
                                    else:
                                        vcf_not_c += 1
                                #print ch, pct[0],pct[1],pct[2]
            acc_pct.extend([it[1] for it in acc_pct_tr])
    print cb_corr_cnt.most_common(30)
    print vcf_c, vcf_not_c
    cov_bin_out = open(cov_bin_f,"w")
    for cbs in cb_corr_cnt:
        cov_bin_out.write("{},{},{}\n".format(cbs[0],cbs[1],cb_corr_cnt[cbs]))
    #pct_bin, pt = np.histogram(acc_pct, bins=500, range=(0, 1))
    #cov_bin_out = open(cov_bin_f,"w")
    #for ix in range(500):
    #    cov_bin_out.write("{},{}\n".format(pt[ix],pct_bin[ix]))
    #print pct_bin
    #print pt


def bam_allele_coverage(bam_in, chr_to_blocks, fa_f, cov_bin_f, vcf_f, cb_seq_dict, min_cnt=100,min_cov=50 ):
    c2i = {"A":0, "C":1, "G":2, "T":3}  # four array.arrays of the same length in order A C G T
    fa_dict={}
    vcf_dict={}
    acc_pct = []
    cb_seq_set = set(cb_seq_dict.keys())
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    #vcf_in = pysam.VariantFile(vcf_f)
    cb_corr_cnt = Counter()
    for ch in chr_to_blocks:
        print ch
        if ch != "chr15":
            continue
        homo_dict = find_homo_regions(fa_dict[ch], chr_to_blocks[ch])
        for ith, bl in enumerate(chr_to_blocks[ch]):
            cnt = bamfile.count(ch, bl.s, bl.e)
            try:
                vcf_dict = dict((it.pos-1, it) for it in vcf_in.fetch(ch[3:], bl.s, bl.e))
            except:
                print ch[3:], "not in vcf.  ",ch
            if cnt < min_cnt:
                continue
            cov = bamfile.count_coverage(ch, bl.s, bl.e,
            quality_threshold=0)  # four array.arrays of the same length in order A C G T
            for v_pos in vcf_dict:
                if v_pos-bl.s>= len(cov[0]):
                    print "SNP position exceed limit.",v_pos-bl.s,len(cov[0])
                    continue
                freq = (cov[0][v_pos-bl.s],cov[1][v_pos-bl.s],cov[2][v_pos-bl.s],cov[3][v_pos-bl.s])
                if sum(freq)<min_cov:
                    continue
                tmp_atcg_set = {}
                for pileupcolumn in bamfile.pileup(ch, v_pos, v_pos+1,truncate=True, min_base_quality=0,ignore_overlaps=False,max_depth=20000):
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            cb_seq, umi_seq = pileupread.alignment.query_name.split("#")[0].split("_")
                            if cb_seq in cb_seq_set:
                                tmp_atcg_set.setdefault(pileupread.alignment.query_sequence[pileupread.query_position],set()).add(cb_seq)
                bs = tmp_atcg_set.keys()
                for ab in range(len(bs)-1):
                    for ab1 in range(ab+1,len(bs)):
                        tmp_set = tmp_atcg_set[bs[ab]] - tmp_atcg_set[snd_b]  # x not y
                        if len(tmp_set)>1:
                            update_corr_cnt(list(tmp_set), cb_corr_cnt)
                        tmp_set = tmp_atcg_set[snd_b] - tmp_atcg_set[bs[ab]]  # y not x
                        if len(tmp_set)>1:
                            update_corr_cnt(list(tmp_set), cb_corr_cnt)
    print cb_corr_cnt.most_common(30)
    cov_bin_out = open(cov_bin_f,"w")
    for cbs in cb_corr_cnt:
        cov_bin_out.write("{},{},{}\n".format(cbs[0],cbs[1],cb_corr_cnt[cbs]))



def get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, out_dir, cb_seq_dict, bam_short, known_position_dict, min_cov=100, report_pct=(0.15,0.85)):
    c2i = {"A":0, "C":1, "G":2, "T":3}  # four array.arrays of the same length in order A C G T
    fa_dict={}
    acc_pct = []
    REF_cnt_dict = {}
    ALT_cnt_dict = {}
    cb_seq_set = set(cb_seq_dict.keys())
    reporting_summary = []
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    if bam_short is not None:
        bam_s = pysam.AlignmentFile(bam_short, "rb")
    cb_corr_cnt = Counter()
    for ch in chr_to_blocks:
        print ch
        homo_dict = find_homo_regions(fa_dict[ch], chr_to_blocks[ch])
        for ith, bl in enumerate(chr_to_blocks[ch]):
            tmp_bl_flat = get_gene_flat({"NNN":bl.transcript_list}, transcript_to_exon)
            for ex in tmp_bl_flat["NNN"]:
                cnt = bamfile.count(ch, ex[0], ex[1])
                if cnt < min_cov:
                    continue
                cov = bamfile.count_coverage(ch, ex[0], ex[1],
                quality_threshold=0)  # four array.arrays of the same length in order A C G T
                if len(cov[0])<20:
                    continue  # ignore tiny exons
                for i in range(5, len(cov[0])-5):  # ignore the bases at the beginning and the end (close to splicing site)
                    tot =  float(cov[0][i]+cov[1][i]+cov[2][i]+cov[3][i])
                    v_pos = ex[0]+i
                    if tot>min_cov and (fa_dict[ch][v_pos]!="N"):
                        freq = cov[c2i[fa_dict[ch][v_pos]]][i]/tot
                        acc_pct.append(freq)
                        base_freq = [("A",cov[0][i]),("C",cov[1][i]),("G",cov[2][i]),("T",cov[3][i])]
                        base_freq.sort(key=lambda x:x[1],reverse=True)
                        if v_pos == 63318364:
                            print base_freq
                        ALT = [it[0] for it in base_freq if it[0] != fa_dict[ch][v_pos]][0] # the most enriched ALT allele
                        alt_freq = cov[c2i[ALT]][i]/tot
                        if (report_pct[0]< alt_freq < report_pct[1]) or ((ch,v_pos) in known_position_dict):
                            tmp_atcg_set = {}
                            if bam_short is not None:
                                try:
                                    cov_s = bam_s.count_coverage(ch, v_pos, v_pos+1, quality_threshold=20)
                                    s_tot = cov_s[0][0]+cov_s[1][0]+cov_s[2][0]+cov_s[3][0]
                                    if s_tot> (min_cov/2):
                                        s_freq = cov_s[c2i[fa_dict[ch][v_pos]]][0]/float(s_tot)
                                    else:
                                        s_freq = -1
                                except:
                                    s_freq = -1
                            else:
                                s_freq = -1
                            seq_ent = seq_entropy(fa_dict[ch][(v_pos-10):(v_pos+10)])
                            indel_freq = -1
                            if ((ch,v_pos) in known_position_dict) or ((ex[0]+i not in homo_dict) and (seq_ent > 1) and (s_freq==-1 or (0.05<s_freq<0.95))):
                                for pileupcolumn in bamfile.pileup(ch, v_pos, v_pos+1,truncate=True, min_base_quality=0,ignore_overlaps=False,max_depth=20000):
                                    c_keep = 0
                                    c_del = 0
                                    for pileupread in pileupcolumn.pileups:
                                        if not pileupread.is_del:
                                            if not pileupread.is_refskip:
                                                c_keep += 1
                                                cb_seq, umi_seq = pileupread.alignment.query_name.split("#")[0].split("_")
                                                if cb_seq in cb_seq_set:
                                                    tmp_atcg_set.setdefault(pileupread.alignment.query_sequence[pileupread.query_position],Counter())[cb_seq] += 1
                                                    #tmp_set[cb_seq] += 1
                                                    if pileupread.alignment.query_sequence[pileupread.query_position] == fa_dict[ch][v_pos]:
                                                        REF_cnt_dict.setdefault((ch, v_pos),[]).append(cb_seq)
                                                    if pileupread.alignment.query_sequence[pileupread.query_position] == ALT:
                                                        ALT_cnt_dict.setdefault((ch, v_pos),[]).append(cb_seq)
                                        else:
                                            if not pileupread.is_refskip:
                                                c_del += 1
                                indel_freq = c_del/float(c_keep+c_del)
                                tmp_set = set()
                                for b in tmp_atcg_set:
                                    tmp_atcg_set[b] = set(it for it in tmp_atcg_set[b] if tmp_atcg_set[b][it]<=2)
                                if (base_freq[0][0] in tmp_atcg_set) and (base_freq[1][0] in tmp_atcg_set):
                                    tmp_set.update(tmp_atcg_set[base_freq[0][0]])
                                    tmp_set.update(tmp_atcg_set[base_freq[1][0]])
                                    rv = hypergeom(len(tmp_set), len(tmp_atcg_set[base_freq[0][0]]), len(tmp_atcg_set[base_freq[1][0]]))
                                    hpg_prob = rv.pmf(len(tmp_atcg_set[base_freq[0][0]].intersection(tmp_atcg_set[base_freq[1][0]])))
                                else:
                                    hpg_prob = 1
                                reporting_summary.append((ch, v_pos, fa_dict[ch][v_pos], ALT, freq, s_freq, hpg_prob, seq_ent, indel_freq))
    print "number:", len(reporting_summary)
    subfolder_name = "mutation"
    if not os.path.exists(os.path.join(out_dir,subfolder_name)):
        os.makedirs(os.path.join(out_dir,subfolder_name))
    with gzip.open(os.path.join(out_dir,subfolder_name,"ref_cnt.csv.gz"),"wb") as ref_cnt_f:
        ref_cnt_f.write("chr,position,"+",".join(cb_seq_dict.keys())+"\n")  # write header
        for p in REF_cnt_dict:
            tmp_c = Counter(REF_cnt_dict[p])
            ref_cnt_f.write("{},{},".format(p[0],p[1])+",".join( str(tmp_c[it]) for it in cb_seq_dict.keys() )+"\n" )
    with gzip.open(os.path.join(out_dir,subfolder_name,"alt_cnt.csv.gz"),"wb") as alt_cnt_f:
        alt_cnt_f.write("chr,position,"+",".join(cb_seq_dict.keys())+"\n")  # write header
        for p in ALT_cnt_dict:
            tmp_c = Counter(ALT_cnt_dict[p])
            alt_cnt_f.write("{},{},".format(p[0],p[1])+",".join( str(tmp_c[it]) for it in cb_seq_dict.keys() )+"\n" )
    with gzip.open(os.path.join(out_dir,subfolder_name,"allele_stat.csv.gz"),"wb") as al_stat:
        al_stat.write("chr,position,REF,ALT,REF_frequency,REF_frequency_in_short_reads,hypergeom_test_p_value,sequence_entrophy,INDEL_frequency\n")  # write header
        for rec in reporting_summary:
            al_stat.write(",".join( str(it) for it in rec )+"\n" )
    pct_bin, pt = np.histogram(acc_pct, bins=500, range=(0, 1))
    with open(os.path.join(out_dir,subfolder_name,"freq_summary.csv"),"w") as cov_bin_out:
        for ix in range(500):
            cov_bin_out.write("{},{}\n".format(pt[ix],pct_bin[ix]))


def get_mito_SNV_table(bam_in, fa_f, out_dir, cb_seq_dict, bam_short, ch="chrM", min_cov=1000, report_pct=(0.15,0.85)):
    c2i = {"A":0, "C":1, "G":2, "T":3}  # four array.arrays of the same length in order A C G T
    fa_dict={}
    acc_pct = []
    REF_cnt_dict = {}
    ALT_cnt_dict = {}
    cb_seq_set = set(cb_seq_dict.keys())
    reporting_summary = []
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    bl = namedtuple("bl", ["s","e"])
    tmp_bl = bl(1, len(fa_dict[ch])-1)
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    bam_s = pysam.AlignmentFile(bam_short, "rb")
    cb_corr_cnt = Counter()
    homo_dict = find_homo_regions(fa_dict[ch], [tmp_bl])
    cnt = bamfile.count(ch, 0, len(fa_dict[ch]))
    cov = bamfile.count_coverage(ch, 0, len(fa_dict[ch]),
    quality_threshold=0)  # four array.arrays of the same length in order A C G T
    for i in range(5, len(cov[0])-5):  # ignore the bases at the beginning and the end
        tot =  float(cov[0][i]+cov[1][i]+cov[2][i]+cov[3][i])
        if tot>min_cov and fa_dict[ch][i] != "N":
            v_pos = i
            freq = cov[c2i[fa_dict[ch][v_pos]]][i]/tot
            acc_pct.append(freq)
            base_freq = [("A",cov[0][i]),("C",cov[1][i]),("G",cov[2][i]),("T",cov[3][i])]
            base_freq.sort(key=lambda x:x[1],reverse=True)
            ALT = [it[0] for it in base_freq if it[0] != fa_dict[ch][v_pos]][0] # the most enriched ALT allele
            alt_freq = cov[c2i[ALT]][i]/tot
            if report_pct[0]< alt_freq < report_pct[1]:
                tmp_atcg_set = {}
                try:
                    cov_s = bam_s.count_coverage(ch, v_pos, v_pos+1, quality_threshold=20)
                    s_tot = cov_s[0][0]+cov_s[1][0]+cov_s[2][0]+cov_s[3][0]
                    if s_tot> (min_cov):
                        s_freq = cov_s[c2i[fa_dict[ch][v_pos]]][0]/float(s_tot)
                    else:
                        s_freq = -1
                except:
                    s_freq = -1
                seq_ent = seq_entropy(fa_dict[ch][(v_pos-10):(v_pos+10)])
                indel_freq = -1
                if (0+i not in homo_dict) and (seq_ent > 1) and (s_freq==-1 or (0.05<s_freq<0.95)):
                    for pileupcolumn in bamfile.pileup(ch, v_pos, v_pos+1,truncate=True, min_base_quality=0,ignore_overlaps=False,max_depth=1000000):
                        mean_base_q = pileupcolumn.get_mapping_qualities()
                        mean_base_q = sum(mean_base_q)/float(len(mean_base_q))
                        c_keep = 0
                        c_del = 0
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del:
                                if not pileupread.is_refskip:
                                    c_keep += 1
                                    cb_seq, umi_seq = pileupread.alignment.query_name.split("#")[0].split("_")
                                    if cb_seq in cb_seq_set:
                                        tmp_atcg_set.setdefault(pileupread.alignment.query_sequence[pileupread.query_position],Counter())[cb_seq] += 1
                                        #tmp_set[cb_seq] += 1
                                        if pileupread.alignment.query_sequence[pileupread.query_position] == fa_dict[ch][v_pos]:
                                            REF_cnt_dict.setdefault((ch, v_pos),[]).append(cb_seq)
                                        if pileupread.alignment.query_sequence[pileupread.query_position] == ALT:
                                            ALT_cnt_dict.setdefault((ch, v_pos),[]).append(cb_seq)
                            else:
                                if not pileupread.is_refskip:
                                    c_del += 1
                    indel_freq = c_del/float(c_keep+c_del)
                    tmp_set = set()
                    for b in tmp_atcg_set:
                        tmp_atcg_set[b] = set(it for it in tmp_atcg_set[b] if tmp_atcg_set[b][it]<=2)
                    if (base_freq[0][0] in tmp_atcg_set) and (base_freq[1][0] in tmp_atcg_set):
                        tmp_set.update(tmp_atcg_set[base_freq[0][0]])
                        tmp_set.update(tmp_atcg_set[base_freq[1][0]])
                        rv = hypergeom(len(tmp_set), len(tmp_atcg_set[base_freq[0][0]]), len(tmp_atcg_set[base_freq[1][0]]))
                        hpg_prob = rv.pmf(len(tmp_atcg_set[base_freq[0][0]].intersection(tmp_atcg_set[base_freq[1][0]])))
                    else:
                        hpg_prob = 1
                    reporting_summary.append((ch, v_pos, fa_dict[ch][v_pos], ALT, freq, s_freq, hpg_prob, seq_ent, indel_freq, mean_base_q))
    if not os.path.exists(os.path.join(out_dir,"mutation")):
        os.makedirs(os.path.join(out_dir,"mutation"))
    with open(os.path.join(out_dir,"mutation","MT_ref_cnt.csv"),"w") as ref_cnt_f:
        ref_cnt_f.write("chr,position,"+",".join(cb_seq_dict.keys())+"\n")  # write header
        for p in REF_cnt_dict:
            tmp_c = Counter(REF_cnt_dict[p])
            ref_cnt_f.write("{},{},".format(p[0],p[1])+",".join( str(tmp_c[it]) for it in cb_seq_dict.keys() )+"\n" )
    with open(os.path.join(out_dir,"mutation","MT_alt_cnt.csv"),"w") as alt_cnt_f:
        alt_cnt_f.write("chr,position,"+",".join(cb_seq_dict.keys())+"\n")  # write header
        for p in ALT_cnt_dict:
            tmp_c = Counter(ALT_cnt_dict[p])
            alt_cnt_f.write("{},{},".format(p[0],p[1])+",".join( str(tmp_c[it]) for it in cb_seq_dict.keys() )+"\n" )
    with open(os.path.join(out_dir,"mutation","MT_allele_stat.csv"),"w") as al_stat:
        al_stat.write("chr,position,REF,ALT,REF_frequency,REF_frequency_in_short_reads,hypergeom_test_p_value,sequence_entrophy,INDEL_frequency,mean_base_quality\n")  # write header
        for rec in reporting_summary:
            al_stat.write(",".join( str(it) for it in rec )+"\n" )
    pct_bin, pt = np.histogram(acc_pct, bins=500, range=(0, 1))
    with open(os.path.join(out_dir,"mutation","MT_freq_summary.csv"),"w") as cov_bin_out:
        for ix in range(500):
            cov_bin_out.write("{},{}\n".format(pt[ix],pct_bin[ix]))

if __name__ == "__main__":
    known_position_dict = {("chr18",63318364):0}
    fa_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCh38.primary_assembly.genome.fa"
    gff_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.v33.annotation.gff3"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    # CLL141 capture
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/Thijssen_count80/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL141-CLL-cells_S8_Rsubread.sorted.bam"
    iso_dir = "/stornext/Genomics/data/CLL_venetoclax/single_cell_data/capture_test/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    """
    ### CLL141
    cb_seq_dict = dict( (it.strip().split(",")[1], it.strip().split(",")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Rachel_scRNA_Aug19/cluster_barcode_anno_lib20.csv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL141-CLL-cells_S8_Rsubread.sorted.bam"
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Rachel_scRNA_Aug19/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    
    ### CLL267
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/patient2/Thijssen2_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short="/stornext/General/data/user_managed/grpu_mritchie_1/hongkePeng/Rachel/all_fastq/CLL267_S4_Rsubread.sorted.bam"
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL267/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    
    ### CLL318
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/CLL318/CLL318_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL318"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    ### CLL152
    
    cb_seq_dict = dict( (it.strip().split(",")[1], it.strip().split(",")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/CLL152/CLL152_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL152/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    
    ### CLL153
    cb_seq_dict = dict( (it.strip().split(",")[0], it.strip().split(",")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Illumina_data/CLL153/cellranger_code/CLL153_count20/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL153/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)


    ### scmix1
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/scbench_5cellline_10x/10percent_cellranger/filtered_gene_bc_matrices/hg38/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_out_8"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)

    ### scmix2
    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/SCmixologyV3/luyiT_10X_260319/cellmix_Lib10/outs/filtered_feature_bc_matrix/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION_April19/v3_long_read/isoform_all"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    """

    """
    #mouse
    known_position_dict = {}
    fa_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/GRCm38.primary_assembly.genome.fa"
    gff_f="/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.vM24.annotation.gff3"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)

    cb_seq_dict = dict( (it.strip().split("-")[0], it.strip().split("-")[0]) for it in open("/stornext/General/data/user_managed/grpu_mritchie_1/JamesRyall/10X/AGRF_CAGRF18671_CD27KANXX_cellranger/MuSC_10P_cellranger/filtered_gene_bc_matrices/mm10/barcodes.tsv"))
    bam_short=None
    iso_dir = "/stornext/General/data/user_managed/grpu_mritchie_1/JamesRyall/PromethION/isoform_out"
    bam_in=os.path.join(iso_dir, "align2genome.bam")
    #gff_f=os.path.join(iso_dir, "isoform_annotated.gff3")
    #chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    #gene_dict = get_gene_flat(gene_to_transcript,transcript_to_exon)
    #chr_to_blocks = get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript)
    get_all_SNV_table(bam_in, chr_to_blocks, transcript_to_exon, fa_f, iso_dir, cb_seq_dict, bam_short,known_position_dict)
    """