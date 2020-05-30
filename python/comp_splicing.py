from itertools import izip
from bisect import bisect_right
from collections import Counter
from parse_gene_anno import parse_gff_tree
import os

# https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

def iv_overlap(iv1, iv2):
    return max(0, min(iv1[1],iv2[1])- max(iv2[0],iv1[0]))


def exon_overlap(exons1, exons2):
    total = 0
    for e in pairwise(exons1):
        total += sum(iv_overlap(e,e2) for e2 in pairwise(exons2))
    return total


def if_exon_contains(s1, s2, max_tol):
    """
    if s2 in s1
    searching for exact match
    """
    if len(s2)==2:  # ignore single exon transcripts
        return False
    if s2[1] not in s1:  # first splicing site
        return False
    fs = s1.index(s2[1])
    if fs == 0 or (s2[0]-s1[fs-1]) < -max_tol:
        return False  # left not within s1
    for i in range(2, len(s2)-1):
        if fs+i-1 > len(s1)-1:
            return False
        if s1[fs+i-1] != s2[i]:
            return False
    if fs+len(s2)-2 > len(s1)-1:
        return False
    if (s2[-1]-s1[fs+len(s2)-2]) > max_tol:
        return False  # right beyond s1
    return True


def make_CAGE_dict(cage_f):
    cage_dict = {}
    for l in open(cage_f):
        items = l.split("\t")
        ch = items[0]
        st = int(items[1])
        en = int(items[2])
        cage_dict.setdefault(ch,[]).append((st, en))
    for ch in cage_dict:
        cage_dict[ch].sort(key=lambda x: x[0])
        new_iv = [cage_dict[ch][0]]
        for i in range(1,len(cage_dict[ch])):
            if cage_dict[ch][i][0]>new_iv[-1][1]:
                new_iv.append(cage_dict[ch][i])
            else:
                new_iv[-1] = (new_iv[-1][0], max(new_iv[-1][1], cage_dict[ch][i][1]) )
        cage_dict[ch] = new_iv
        if not all(cage_dict[ch][i-1][1]<cage_dict[ch][i][0] for i in range(1,len(cage_dict[ch]))):
            print ch,"overlapping exists in bed file"
    return cage_dict


def get_cage_coverage(isoform_gff, cage_f):
    cage_dict = make_CAGE_dict(cage_f)
    cage_cov_dict = {}
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(isoform_gff)
    for ch in chr_to_gene:
        if ch not in cage_dict:
            print(ch,"not in CAGE annotation file.")
            continue
        cage_cov_dict[ch] = []
        cage_tmp = dict((it[0],it) for it in cage_dict[ch])
        cage_left = [it[0] for it in cage_dict[ch]]
        for ge in chr_to_gene[ch]:
            for tr in gene_to_transcript[ge]:
                tss_pos = transcript_dict[tr].start if transcript_dict[tr].strand=="+" else transcript_dict[tr].end
                min_idx = max(0,bisect_right(cage_left, tss_pos)-1)
                #min_idx = min(range(len(cage_dict[ch])), key=lambda i: min(abs(cage_dict[ch][i][0]-tss_pos),abs(cage_dict[ch][i][1]-tss_pos)) )
                if tss_pos>=cage_dict[ch][min_idx][0] and tss_pos<=cage_dict[ch][min_idx][1]:
                    cage_cov_dict[ch].append((tr,tss_pos,0))
                elif tss_pos<cage_dict[ch][min_idx][0]:
                    if min_idx==0:
                        cage_cov_dict[ch].append((tr,tss_pos,cage_dict[ch][min_idx][0]-tss_pos))
                    else:
                        cage_cov_dict[ch].append((tr,tss_pos,min(tss_pos-cage_dict[ch][min_idx-1][1], cage_dict[ch][min_idx][0]-tss_pos ) ))
                else:  # tss_pos>cage_dict[ch][min_idx][1]
                    if min_idx==len(cage_dict[ch])-1:
                        cage_cov_dict[ch].append((tr,tss_pos,tss_pos-cage_dict[ch][min_idx][1]))
                    else:
                        cage_cov_dict[ch].append((tr,tss_pos,min(tss_pos-cage_dict[ch][min_idx][1],cage_dict[ch][min_idx+1][0]-tss_pos )))
    return cage_cov_dict


def comp_two_splicing(s1, s2):
    """
    compare s2 to s1.
    no_tr_overlap: 
        transcript not overlap
    no_splice_overlap:
        transcript overlap but no common splice site
    """
    typ = []
    if s1[0]> s2[-1] or s2[0] > s1[-1]:
        typ.append("no_tr_overlap")
        return typ
    elif all(x not in s1 for x in s2[1:-1]):
        typ.append("no_splice_overlap")
        return typ
    if len(s1)==len(s2) and all(x1 == x2 for x1,x2 in zip(s1[1:-1],s2[1:-1])):
        typ.append("full_splice_match")
        return typ
    if all(x in s1 for x in s2[1:-1]): # all splicing site in s2 can be found in s1
        s_idx_list = [s1.index(x) for x in s2[1:-1]]
        if all(s_idx_list[x]-s_idx_list[x-1]==1 for x in range(1,len(s_idx_list))):  # no exon skipping
            if (s_idx_list[0]-2>0 and s2[0]>s1[s_idx_list[0]-2]) or s_idx_list[0]==1:
                if (s_idx_list[-1]+2<len(s1) and s2[-1]<s1[s_idx_list[-1]+2]) or s_idx_list[-1]==len(s1)-2:
                    typ.append("incomplete_splice_match")
                else:
                    typ.append("intron_retention")
            else:
                typ.append("intron_retention")
        else:
            tmp_overlap = [[iv_overlap(e,e2)>0 for e2 in pairwise(s1)] for e in pairwise(s2)]
            if any(sum(it)>1 for it in tmp_overlap):  # exon in s2 overlap to more than one exon in s1
                typ.append("intron_retention")
            else:
                typ.append("skip_exon")
    else:
        typ.append("not_all_splice_match") # some splicing site in s2 cannot be found in s1
        s_idx_list = [s1.index(x) if x in s1 else -1 for x in s2[1:-1]]
        if s_idx_list.count(-1)==1:
            typ.append("diff_one_splice_site")
        if s_idx_list[0]==-1 and all(x>=0 for x in s_idx_list[1:]):
            typ.append("alt_first_exon")
        elif s_idx_list[-1]==-1 and all(x>=0 for x in s_idx_list[:-1]):
            typ.append("alt_last_exon")  # all other splice site are the same except last splice site
        if s1[1]==s2[1]:
            typ.append("same_fist_splice")
        if s1[-2]==s2[-2]:
            typ.append("same_last_splice")
        if ("same_fist_splice" in typ) and ("same_last_splice" in typ):
            if all(x in s2 for x in s1[1:-1]): # all splicing site in s1 can be found in s2
                typ.append("exon_inclusion")
            else:
                s1_pair = [e for e in pairwise(s1)]
                s2_pair = [e for e in pairwise(s2)]
                if len(s1_pair)>2 and len(s2_pair)>2:
                    for s2_p in s2_pair[1:-1]:
                        if s2_p not in s1_pair and any(iv_overlap(e1,s2_p)>0 for e1 in s1_pair):  # partialy overlap
                            typ.append("exon_overlap_diff_splice_site")
                            break 
                    if "exon_overlap_diff_splice_site" not in typ:
                        typ.append("mutually_exclusive_exon")
                else:
                    typ.append("others")
        if s1[1]<s2[0] or s2[1]<s1[0]:  # no overlap in first exon
            typ.append("no_first_exon_overlap")
        if s1[-2]>s2[-1] or s2[-2]>s1[-1]:  # no overlap in last exon
            typ.append("no_last_exon_overlap")
        if len(typ)==1:
            typ.append("todo")
        #for e in pairwise(s2)
    return typ


def get_splice_expr(isoform_gff, fsm_annotation, res_csv, tr_kept=10):
    typ_cnt = Counter()
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(isoform_gff)
    transcript_to_splice = {}
    for tr in transcript_to_exon:  # convert to list
        transcript_to_splice[tr] = []
        for ex in transcript_to_exon[tr]:
            transcript_to_splice[tr].append(ex[0])
            transcript_to_splice[tr].append(ex[1])
    tr_dict = {}
    fsm2tr = {}
    fsm_cnt = Counter()
    gene2fsm = {}
    for line in open(fsm_annotation):
        if "gene_id" in line:
            continue
        its = line.strip().split(",")  # 0:tr_id, 1:gene_id, 2:FSM_id, 3:in_ref, 4:cnt
        if its[0] not in transcript_to_exon:  # filtered out in gff (low abundance)
            continue
        tr_dict[its[0]] = (its[1],its[2],int(its[4]))
        fsm_cnt[its[2]] += int(its[4])
        fsm2tr.setdefault(its[2],[]).append(its[0])
        gene2fsm.setdefault(its[1],[]).append(its[2])
    for ge in gene2fsm:
        gene2fsm[ge] = [(it, fsm_cnt[it]) for it in list(set(gene2fsm[ge]))]
        if len(gene2fsm[ge])==1:
            continue
        gene2fsm[ge].sort(key=lambda x: x[1],reverse=True)
        if len(gene2fsm[ge])>tr_kept:
            gene2fsm[ge] = gene2fsm[ge][:tr_kept]
    out_f = open(res_csv,"w")
    for ge in gene2fsm:
        if len(gene2fsm[ge])==1:
            continue
        for i in range(len(gene2fsm[ge])):
            for j in range(len(gene2fsm[ge])):
                if i==j:
                    continue
                if gene2fsm[ge][i][0] in transcript_to_splice:
                    tr1 = gene2fsm[ge][i][0]
                else:
                    tr1 = fsm2tr[gene2fsm[ge][i][0]][0]
                if gene2fsm[ge][j][0] in transcript_to_splice:
                    tr2 = gene2fsm[ge][j][0]
                else:
                    tr2 = fsm2tr[gene2fsm[ge][j][0]][0]
                typ = comp_two_splicing(transcript_to_splice[tr1], transcript_to_splice[tr2])
                typ = "::".join(typ)
                typ_cnt[typ] += 1
                out_f.write("{},{},{}\n".format(gene2fsm[ge][i][0],gene2fsm[ge][j][0], typ))
    out_f.close()
    for i in typ_cnt.most_common(5):
        print "\t",i



if __name__ == "__main__":
    data_dir_list = ["/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoform_all",
    "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION_April19/v3_long_read/isoform_all",
    "/stornext/General/data/user_managed/grpu_mritchie_1/JamesRyall/PromethION/isoform_out",
    "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/Rachel_scRNA_Aug19/isoform_out",
    "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL267/isoform_out",
    "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL152/isoform_out",
    "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL153/isoform_out",
    "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL318/isoform_out"]
    for data_dir in data_dir_list:
        print data_dir
        isoform_gff = os.path.join(data_dir,"isoform_annotated.filtered.gff3")
        fsm_annotation = os.path.join(data_dir,"isoform_FSM_annotation.csv")
        res_csv = os.path.join(data_dir,"fsm_splice_comp.csv")
        get_splice_expr(isoform_gff, fsm_annotation, res_csv)
    """
    cage_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/script/SQANTI2_scripts/hg38.cage_peak_phase1and2combined_coord.bed"
    isoform_gff = "/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL152/isoform_out/isoform_annotated.gff3"
    ref_gff = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/gencode.v29.annotation.gff3"
    cage_cov_dict = get_cage_coverage(isoform_gff, cage_f)
    cage_cov_dict_ref = get_cage_coverage(ref_gff, cage_f)
    with open("/stornext/General/data/user_managed/grpu_mritchie_1/RachelThijssen/sclr_data/CLL152/isoform_out/cage_coverage.csv","w") as out_f:
        out_f.write("transcript_id,chromosome,TSS,dis_to_nearest_cage_peak\n")
        for ch in cage_cov_dict:
            for i in range(len(cage_cov_dict[ch])):
                out_f.write("{},{},{},{}\n".format(cage_cov_dict[ch][i][0],ch,cage_cov_dict[ch][i][1],cage_cov_dict[ch][i][2]))
            for i in range(len(cage_cov_dict_ref[ch])):
                out_f.write("{},{},{},{}\n".format(cage_cov_dict_ref[ch][i][0],ch,cage_cov_dict_ref[ch][i][1],cage_cov_dict_ref[ch][i][2]))
    """