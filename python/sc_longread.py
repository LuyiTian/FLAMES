# find isoforms in longread data
import pysam
#from BCBio import GFF # never stop running. is there a bug?

####### from: https://techoverflow.net/2013/11/30/a-simple-gff3-parser-in-python/
from collections import namedtuple, Counter
import gzip
from urlparse import urlparse
import os
import numpy as np
import bisect
import random
import copy
from itertools import izip
from parse_gene_anno import parse_gff_tree


def take_closest(l, num):
    return min(l, key=lambda x: abs(x-num))


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
        total += sum(iv_overlap(e,e2) for e2 in exons2)
    return total


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


## pysam cigar specification:

# M   BAM_CMATCH  0
# I   BAM_CINS    1
# D   BAM_CDEL    2
# N   BAM_CREF_SKIP   3
# S   BAM_CSOFT_CLIP  4
# H   BAM_CHARD_CLIP  5
# P   BAM_CPAD    6
# =   BAM_CEQUAL  7
# X   BAM_CDIFF   8
# B   BAM_CBACK   9


def smooth_cigar(cigar_tup,thr=10):
    new_cigar = [list(cigar_tup[0])]
    for ix in range(1,len(cigar_tup)):
        if new_cigar[-1][0] != 0:
            new_cigar.append(list(cigar_tup[ix]))
        elif cigar_tup[ix][0] == 0:  # merge matched reads
            new_cigar[-1][1] = new_cigar[-1][1] + cigar_tup[ix][1] 
        elif cigar_tup[ix][0] == 2:
            if cigar_tup[ix][1]<=thr:
                new_cigar[-1][1] = new_cigar[-1][1] + cigar_tup[ix][1] # add deletion to match
            else:
                new_cigar.append(list(cigar_tup[ix]))  # big deletion
        elif cigar_tup[ix][0] == 1:
            if cigar_tup[ix][1]>thr:
                new_cigar.append(list(cigar_tup[ix]))  # big insertion
            else:
                continue
        else:
            new_cigar.append(list(cigar_tup[ix]))
    return new_cigar


def get_all_site(transcript_to_junctions, tr_list):
    all_site = set()
    for t in tr_list:
        all_site.add(transcript_to_junctions[t]["left"])
        all_site.add(transcript_to_junctions[t]["right"])
        all_site.update(transcript_to_junctions[t]["junctions"])
    all_site = list(all_site)
    all_site.sort()
    return all_site


def get_splice_site(transcript_to_junctions, tr_list):
    all_site = set()
    for t in tr_list:
        all_site.update(transcript_to_junctions[t]["junctions"])
    all_site = list(all_site)
    all_site.sort()
    return all_site


def get_TSS_TES_site(transcript_to_junctions, tr_list):
    all_site = {"left":[],"right":[]}
    for t in tr_list:
        if len(all_site["left"]) > 0:
            if abs(take_closest(all_site["left"], transcript_to_junctions[t]["left"])-transcript_to_junctions[t]["left"]) > 5:
                all_site["left"].append(transcript_to_junctions[t]["left"])
        else:
            all_site["left"].append(transcript_to_junctions[t]["left"])
        if len(all_site["right"]) > 0:
            if abs(take_closest(all_site["right"], transcript_to_junctions[t]["right"])-transcript_to_junctions[t]["right"]) > 5:
                all_site["right"].append(transcript_to_junctions[t]["right"])
        else:
            all_site["right"].append(transcript_to_junctions[t]["right"])
    all_site["left"].sort()
    all_site["right"].sort()
    return all_site


def get_exon_sim_pct(exons1, exons2):
    """
    get percentage of coverage between two transcripts
    """
    def pos_overlap(pos1, pos2):
        if pos1[1]<=pos2[0] or pos1[0]>=pos2[1]:
            return 0
        else:
            return min(pos1[1],pos2[1]) - max(pos1[0],pos2[0])
    e1_len = sum(p[1]-p[0] for p in pairwise(exons1))
    e2_len = sum(p[1]-p[0] for p in pairwise(exons2))
    total = 0
    for e in pairwise(exons1):
        total += sum(pos_overlap(e,e2) for e2 in pairwise(exons2))
    return float(total)/max(e1_len, e2_len)


def is_exon_similar(ex1, ex2, thr):
    if len(ex1) != len(ex2):
        return False
    cum_diff = 0
    for i in range(len(ex1)):
        cum_diff += abs(ex1[i][0]-ex2[i][0])+abs(ex1[i][1]-ex2[i][1])
        if cum_diff > thr:
            return False
    return True


def remove_similar_tr(transcript_dict, gene_to_transcript, transcript_to_exon, thr=10):
    dup_stat = Counter()
    for g in gene_to_transcript:
        if len(gene_to_transcript[g]) < 2:
            continue
        dup_list = []
        for tr_idx in range(len(gene_to_transcript[g])-1):
            for tr2_idx in range(tr_idx+1,len(gene_to_transcript[g])):
                if is_exon_similar(transcript_to_exon[gene_to_transcript[g][tr_idx]],transcript_to_exon[gene_to_transcript[g][tr2_idx]], thr):
                    dup_list.append(tr2_idx)
                    dup_stat["duplicated_transcripts"] += 1
                    #print gene_to_transcript[g][tr2_idx]
        if len(dup_list) > 0:
            dup_list = list(set(dup_list))
            gene_to_transcript[g] = [i for j, i in enumerate(gene_to_transcript[g]) if j not in dup_list]
    print "remove similar transcripts in gene annotation:", dup_stat


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


def find_best_splice_chain(raw_iso, junc_list, max_dist):
    best_match = [-1,0,0]  # index, length-1, starting position
    for ith, junc in enumerate(junc_list):
        i_st = [ix for ix,it in enumerate(junc) if abs(it-raw_iso[1])<max_dist]
        if len(i_st)==1:
            i_st = i_st[0]
            i = i_st+1
            i_end = i_st
            while i<len(junc) and i-i_st+1<len(raw_iso):
                if abs(junc[i]-raw_iso[i-i_st+1])<max_dist:
                    i_end = i
                    i += 1
                else:
                    break
            if i_end-i_st>=best_match[1]:
                best_match[0] = ith
                best_match[1] = i_end-i_st
                best_match[2] = i_st
    if best_match[0]>=0 and best_match[1]>=3:
        updated_iso = list(raw_iso)
        for it in range(best_match[2],best_match[2]+best_match[1]+1):
            updated_iso[it-best_match[2]+1] = junc_list[best_match[0]][it]  # update splicing site
        return updated_iso
    else:
        return list(raw_iso)


def get_gene_flat(gene_to_transcript, transcript_to_exon):
    gene_dict = {}
    for g in gene_to_transcript:
        gene_dict[g] = copy.deepcopy(transcript_to_exon[gene_to_transcript[g][0]])
        if len(gene_to_transcript[g]) > 1:
            for i in range(1, len(gene_to_transcript[g])):
                for j in range(len(transcript_to_exon[gene_to_transcript[g][i]])):
                    j_found = False
                    for k in range(len(gene_dict[g])):
                        if transcript_to_exon[gene_to_transcript[g][i]][j][0]> gene_dict[g][k][1] or \
                        transcript_to_exon[gene_to_transcript[g][i]][j][1]< gene_dict[g][k][0]:
                            continue
                        j_found = True
                        if transcript_to_exon[gene_to_transcript[g][i]][j][1] > gene_dict[g][k][1]:
                            gene_dict[g][k][1] = transcript_to_exon[gene_to_transcript[g][i]][j][1]
                        if transcript_to_exon[gene_to_transcript[g][i]][j][0] < gene_dict[g][k][0]:
                            gene_dict[g][k][0] = transcript_to_exon[gene_to_transcript[g][i]][j][0]
                    if j_found==False:
                        gene_dict[g].append(transcript_to_exon[gene_to_transcript[g][i]][j]) # no overlapping, add new exon
        gene_dict[g].sort(key=lambda x:x[0]) # sort by left most positions
    return gene_dict

class GeneBlocks(object):
    """docstring for GeneBlocks"""
    def __init__(self, start, end, transcript_list, a_gene):
        self.s = start
        self.e = end
        self.transcript_list = transcript_list
        self.gene_to_tr = {}
        self.gene_to_tr[a_gene] = tuple(transcript_list)
    def add_gene(self, start, end, transcript_list, a_gene):
        self.e = max(self.e, end)
        self.transcript_list.extend(transcript_list)
        self.gene_to_tr[a_gene] = tuple(transcript_list)
    def __str__(self):
        print("GeneBlocks({},{}): {} transcripts, {} genes.".format(self.s, self.e, len(self.transcript_list), len(self.gene_to_tr)))
    def __repr__(self):
        return "GeneBlocks({},{}): {} transcripts, {} genes.".format(self.s, self.e, len(self.transcript_list), len(self.gene_to_tr))


def blocks_to_junctions(blocks):
    junctions = {"left":blocks[0][0], "right":blocks[-1][1], "junctions":[]}
    if len(blocks)>1:
        for i in range(1,len(blocks)):
            junctions["junctions"].append(blocks[i-1][1])
            junctions["junctions"].append(blocks[i][0])
        junctions["junctions"] = tuple(junctions["junctions"])
    else:
        junctions["junctions"] = tuple()
    return junctions


def get_gene_blocks(gene_dict, chr_to_gene, gene_to_transcript):
    chr_to_blocks = {}
    for ch in chr_to_gene:
        gene_l = []
        chr_to_blocks[ch] = []
        for g in chr_to_gene[ch]:
            gene_l.append([g, gene_dict[g][0][0], gene_dict[g][-1][1]]) # gene_id, start, end
        gene_l.sort(key=lambda x:x[1]) # sort based on starting position
        chr_to_blocks[ch].append(GeneBlocks(gene_l[0][1],
            gene_l[0][2],
            gene_to_transcript[gene_l[0][0]],
            gene_l[0][0]))
        if len(gene_l) > 1:
            for i in range(1,len(gene_l)):
                if chr_to_blocks[ch][-1].e > gene_l[i][1]:
                    chr_to_blocks[ch][-1].add_gene(gene_l[i][1], gene_l[i][2], gene_to_transcript[gene_l[i][0]], gene_l[i][0])
                else:
                    chr_to_blocks[ch].append(GeneBlocks(gene_l[i][1],
                        gene_l[i][2],
                        gene_to_transcript[gene_l[i][0]],
                        gene_l[i][0]))
    return chr_to_blocks


def parse_bam_intron(bam_in, bam_out, summary_csv, chr_to_blocks, gene_dict):
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    bamfile_o = pysam.AlignmentFile(bam_out, "wb", template=bamfile)
    csv_out = open(summary_csv,"w")
    map_tag="YE"
    gene_tag="GE"
    for ch in chr_to_blocks:
        # print ch
        for bl in chr_to_blocks[ch]:
            it_region = bamfile.fetch(ch, bl.s, bl.e)
            for rec in it_region:
                if rec.is_secondary:
                    continue
                #rec.cigar = smooth_cigar(rec.cigar)
                if len(bl.gene_to_tr)==1:
                    ov_exon = sum(rec.get_overlap(it[0], it[1]) for it in gene_dict[bl.gene_to_tr.keys()[0]])
                    if ov_exon == 0:
                        found = False
                        rec.set_tag(map_tag, 10) # unique mapped only to intron
                        rec.set_tag(gene_tag, bl.gene_to_tr.keys()[0], "Z")
                        bamfile_o.write(rec)
                        csv_out.write("{chr}\t{start}\t{end}\t{name}\t{mapping_stat}\t{strand}\t{gene_name}\n".format(
                            chr=ch,start=rec.reference_start,end=rec.reference_start+rec.reference_length,
                            name=rec.query_name,mapping_stat=10,strand="-" if rec.is_reverse else "+",gene_name=bl.gene_to_tr.keys()[0]))
                    else:
                        found = True
                else:
                    found = False
                    gene_name = ""
                    for g in bl.gene_to_tr:
                        ov_exon = sum(rec.get_overlap(it[0], it[1]) for it in gene_dict[g])
                        ov_gene = rec.get_overlap(gene_dict[g][0][0], gene_dict[g][-1][1])
                        if ov_exon==0 and ov_gene>0: # in genebody, not overlap to exon
                            if not found: # not in exon
                                if gene_name == "": # first time in intron
                                    gene_name = g
                                    mapping_stat = 10
                                else: # not first time in intron, ambiguous
                                    mapping_stat = 11
                        if ov_exon>0:
                            found = True
                    if not found: # not in exon
                        rec.set_tag(map_tag, mapping_stat) # unique(10)/ambiguous(11) mapped only to intron
                        csv_out.write("{chr}\t{start}\t{end}\t{name}\t{mapping_stat}\t{strand}\t{gene_name}\n".format(
                            chr=ch,start=rec.reference_start,end=rec.reference_start+rec.reference_length,
                            name=rec.query_name,mapping_stat=mapping_stat,strand="-" if rec.is_reverse else "+",gene_name=bl.gene_to_tr.keys()[0] if mapping_stat==11 else "."))
                        if mapping_stat==10: # if unique
                            rec.set_tag(gene_tag, gene_name, "Z")
                        bamfile_o.write(rec)
    bamfile.close()
    bamfile_o.close()
    csv_out.close()

class Isoforms(object):
    """docstring for Isoforms"""
    def __init__(self, ch, config):
        self.ch = ch
        self.junction_dict = {}
        self.junction_list = []
        self.lr_pair = {}
        self.left = []
        self.right = []
        self.single_block_dict = {}
        self.single_blocks = []
        self.MAX_DIST=config["MAX_DIST"]
        self.MAX_TS_DIST=config["MAX_TS_DIST"]
        self.MAX_SPLICE_MATCH_DIST=config["MAX_SPLICE_MATCH_DIST"]
        self.Max_site_per_splice = config["Max_site_per_splice"]
        self.Min_sup_cnt=config["Min_sup_cnt"]
        self.min_fl_exon_len = config["min_fl_exon_len"]
        self.Min_sup_pct=config["Min_sup_pct"]
        self.strand_specific = config["strand_specific"]
        self.remove_incomp_reads= config["remove_incomp_reads"] # 0: keep all truncated isoforms. 1: remove truncated reads (modest level). 1: remove truncated reads (unless they really highly expressed)
        self.strand_cnt = {} # 0: reverse strand. 1: forward strand
        self.new_isoforms = {} # new isoforms after matching the annotation
        self.known_isoforms = {} # known isoforms after matching the annotation
        self.raw_isoforms = {} #isoform assembly before matching the annotation
        self.ge_dict = {} # gene id to isoform

    def add_isoform(self, junctions, is_rev):
        """
        junctions: dict{"left":int, "right":int, "junctions":list}
        is_rev: true if mapped to reverse strand
        """
        if len(junctions["junctions"])==0:  # single exon reads
            if len(self.single_block_dict)==0:
                self.__add_one(junctions,is_rev)
            else:
                if (junctions["left"],junctions["right"]) in self.single_block_dict:
                    self.__update_one(junctions,(junctions["left"],junctions["right"]),is_rev)
                else:
                    found = False
                    for j in self.single_block_dict:
                        if abs(j[0]-junctions["left"])<self.MAX_TS_DIST and abs(j[1]-junctions["right"])<self.MAX_TS_DIST:
                            self.__update_one(junctions,j,is_rev)
                            found = True
                            break
                    if not found:
                        self.__add_one(junctions,is_rev)
        else:
            if len(self.junction_dict)==0:
                self.__add_one(junctions,is_rev)
            else:
                if tuple(junctions["junctions"]) in self.junction_dict:
                    self.__update_one(junctions,tuple(junctions["junctions"]),is_rev)
                else:
                    found = False
                    for j in self.junction_dict:
                        if len(j)==len(junctions["junctions"]):
                            if all(abs(it-ij)<self.MAX_DIST for it, ij in zip(j, junctions["junctions"])):
                                self.__update_one(junctions,j,is_rev)
                                found = True
                                break
                    if not found:
                        self.__add_one(junctions,is_rev)

    def __add_one(self, junctions, strand):
        if len(junctions["junctions"])==0:
            self.single_block_dict[tuple((junctions["left"],junctions["right"]))] = len(self.single_blocks)
            self.single_blocks.append([tuple((junctions["left"],junctions["right"]))])
            self.strand_cnt[tuple((junctions["left"],junctions["right"]))]=[]
            self.strand_cnt[tuple((junctions["left"],junctions["right"]))].append(-1*self.strand_specific if strand else 1*self.strand_specific)
        else:
            self.junction_dict[tuple(junctions["junctions"])] = len(self.junction_list)
            self.junction_list.append([tuple(junctions["junctions"])])
            self.left.append(junctions["left"])
            self.right.append(junctions["right"])
            self.lr_pair[tuple(junctions["junctions"])] = []
            self.lr_pair[tuple(junctions["junctions"])].append((junctions["left"], junctions["right"]))
            self.strand_cnt[tuple(junctions["junctions"])]=[]
            self.strand_cnt[tuple(junctions["junctions"])].append(-1*self.strand_specific if strand else 1*self.strand_specific) 
            # by default long read is opposite to the direction of transcript (self.strand_specific=-1)

    def __update_one(self, junctions, st, strand):
        """
        update an existing isoform. st is the key of dict
        """
        if len(junctions["junctions"])==0:
            self.single_blocks[self.single_block_dict[st]].append(
                        tuple((junctions["left"],junctions["right"])))
        else:
            self.junction_list[self.junction_dict[st]].append(tuple(junctions["junctions"]))
            self.left.append(junctions["left"])
            self.right.append(junctions["right"])
            self.lr_pair[st].append((junctions["left"], junctions["right"]))
        self.strand_cnt[st].append(-1*self.strand_specific if strand else 1*self.strand_specific)

    def __len__(self):
        return len(self.junction_dict)+len(self.single_block_dict)

    def update_all_splice(self):
        junction_tmp = {}
        single_block_tmp = {}
        strand_cnt_tmp = {}
        lr_pair_tmp = {}
        for j in self.junction_dict:
            if len(self.junction_list[self.junction_dict[j]]) >= self.Min_sup_cnt:  # only more than `Min_sup_cnt` reads support this splicing
                new_j = tuple((Counter(it[i] for it in self.junction_list[self.junction_dict[j]]).most_common(1)[0][0] for i in range(len(j))))
                junction_tmp[new_j] = self.junction_dict[j]
                strand_cnt_tmp[new_j] =  Counter(self.strand_cnt[j]).most_common(1)[0][0]
                lr_pair_tmp[new_j] = self.lr_pair[j]
        for j in self.single_block_dict:
            if len(self.single_blocks[self.single_block_dict[j]]) >= self.Min_sup_cnt:
                new_j = tuple((Counter(it[i] for it in self.single_blocks[self.single_block_dict[j]]).most_common(1)[0][0] for i in range(len(j))))
                single_block_tmp[new_j] = self.single_block_dict[j]
                strand_cnt_tmp[new_j] =  Counter(self.strand_cnt[j]).most_common(1)[0][0]
        self.single_block_dict = single_block_tmp
        self.junction_dict = junction_tmp
        self.strand_cnt = strand_cnt_tmp
        self.lr_pair = lr_pair_tmp

    def update_start_end(self):
        pass

    def __group_sites(self, l, smooth_window=3, min_thr=9999):
        """
        `smooth_window` should be odd number
        only smooth if more than `min_thr` reads
        """
        l_cnt = Counter(l)
        l_uni = np.array(list(l_cnt.keys()))
        l_uni.sort()
        if len(l)<min_thr:
            return l_uni, l_cnt
        if l_uni[-1]-l_uni[0] < 5: # all in the same place
            return l_uni, l_cnt
        x = np.zeros(l_uni[-1]-l_uni[0]+1+(smooth_window-1)/2) # maybe I should use sparse matrix
        x[l_uni-l_uni[0]] = [l_cnt[it] for it in l_uni]
        w = np.hanning(smooth_window)
        w[(smooth_window-1)/2]=2*w[(smooth_window-1)/2]
        #x = np.array([l_cnt[it] for it in l_uni])
        s=np.r_[x[smooth_window-1:0:-1],x,x[-2:-smooth_window-1:-1]]
        #y=np.convolve(w/w.sum(),s,mode='valid')
        y=np.convolve(w,s,mode='valid')
        for i in range(len(l_uni)):
            l_cnt[l_uni[i]] = y[l_uni[i]-l_uni[0]+(smooth_window-1)/2]
        return l_uni, l_cnt

    def site_stat(self,out_f):
        """
        out_f should be opened file
        """
        bedgraph_fmt = "{_ch}\t{_st}\t{_en}\t{_sc}\n"
        if len(self.left)==0:
            return 0
        l_uni, l_cnt = self.__group_sites(self.left)
        for i in l_uni:
            out_f.write(bedgraph_fmt.format(_ch=self.ch, _st=i, _en=i+1, _sc=l_cnt[i]))
        l_uni, l_cnt = self.__group_sites(self.right)
        for i in l_uni:
            out_f.write(bedgraph_fmt.format(_ch=self.ch, _st=i, _en=i+1, _sc=l_cnt[i]))
        all_splice = []
        for j in self.junction_list:
            for aj in j:
                all_splice.extend([i for i in aj])
        l_uni, l_cnt = self.__group_sites(all_splice)
        for i in l_uni:
            out_f.write(bedgraph_fmt.format(_ch=self.ch, _st=i, _en=i+1, _sc=l_cnt[i]))

    def filter_TSS_TES(self, out_f, known_site=None, fdr_cutoff=0.01):
        bedgraph_fmt = "{_ch}\t{_st}\t{_en}\t{_sc}\n"
        def filter_site(l_cnt):
            mx = np.array(l_cnt.most_common())
            #print mx[:,:5]
            #mx[:,1] = np.log(mx[:,1])
            if mx[0,1]==1:
                return mx[:,]
            elif mx.shape[0]<5:
                return mx[:,]
            trun_num = max(2,int(0.05*mx.shape[0]))
            rate = 1/np.mean(mx[:,1][:(-trun_num)])
            prob = rate*np.exp(-rate*mx[:,1])  # use exponential distribution
            cum_prob = np.cumsum(prob)
            #print trun_num,"~~~~~" , rate, mx[:,1][:3], mx.shape, "~~~", prob[:3]
            if cum_prob[0]>fdr_cutoff: # no p-value smaller the cutoff
                return mx[:,]
            else:
                idx = np.argmax(cum_prob>fdr_cutoff)
                return mx[:idx,]

        def insert_dist(fs, known_site):
            tmp = [fs[0]]
            if len(fs)>1:
                for s in fs[1:]:
                    clo_p = take_closest(tmp, s)
                    if abs(clo_p-s)>self.MAX_TS_DIST/2:
                        tmp.append(s)
            if known_site is not None:
                for s in known_site:
                    clo_p = take_closest(tmp, s)
                    if abs(clo_p-s)>self.MAX_TS_DIST:
                        tmp.append(s)
            return tmp

        if len(Counter(self.left))<self.Min_sup_cnt or len(Counter(self.right))<self.Min_sup_cnt:
            return 0 
        else:
            #left
            fs_l = filter_site(Counter(self.left))
            #fs_l = filter_site(self.__group_sites(self.left)[1])
            if fs_l.shape[0] == 0:
                #print self.ch, "no left:", self.left
                cnt_l = [-99999999]
            else:
                null = [out_f.write(bedgraph_fmt.format(_ch=self.ch, _st=int(fs_l[ix,0]), _en=int(fs_l[ix,0]+1), _sc=fs_l[ix,1])) for ix in range(fs_l.shape[0])]
                cnt_l = insert_dist(fs_l[:,0],known_site["left"] if known_site is not None else None)
                cnt_l.sort()
            fs_r = filter_site(Counter(self.right))
            #fs_r = filter_site(self.__group_sites(self.right)[1])
            if fs_r.shape[0] == 0:
                #print self.ch, "no right:", self.right
                cnt_r = [-99999999]
            else:
                null = [out_f.write(bedgraph_fmt.format(_ch=self.ch, _st=int(fs_r[ix,0]), _en=int(fs_r[ix,0]+1), _sc=fs_r[ix,1])) for ix in range(fs_r.shape[0])]
                cnt_r = insert_dist(fs_r[:,0],known_site["right"] if known_site is not None else None)
                cnt_r.sort()
        for j in self.lr_pair:
            tmp_pair = Counter(self.lr_pair[j]).most_common()
            pair_after_filtering = []
            pair_enrich = []
            for p,_ in tmp_pair:  # first search for common enriched TSS/TES
                cl_p = (take_closest(cnt_l, p[0]), take_closest(cnt_r, p[1]))
                if (abs(cl_p[0]-p[0]) < 0.5*self.MAX_TS_DIST)  and  (abs(cl_p[1]-p[1]) < 0.5*self.MAX_TS_DIST):
                    if len(pair_after_filtering) > 0:
                        dis = abs(take_closest([it[0] for it in pair_after_filtering], cl_p[0]) -cl_p[0])+abs(take_closest([it[1] for it in pair_after_filtering], cl_p[1])-cl_p[1])
                        if dis > self.MAX_TS_DIST:
                            if cl_p[0]<j[0] and cl_p[1]>j[-1]:
                                pair_after_filtering.append(cl_p)
                    else:
                        if cl_p[0]<j[0] and cl_p[1]>j[-1]:
                            pair_after_filtering.append(cl_p)
                    if len(pair_after_filtering) >= self.Max_site_per_splice: # only allow up to `Max_site_per_splice` combinations
                        break
            if len(pair_after_filtering) == 0:  # then search for isoform specific enrichment
                for p, _ in tmp_pair:
                    pair_enrich.append((p, sum(it[1] for it in tmp_pair if abs(it[0][0]-p[0])+abs(it[0][1]-p[1])< self.MAX_TS_DIST)))
                pair_enrich.sort(key=lambda x: x[1], reverse=True)
                for p, _ in pair_enrich:
                    if len(pair_after_filtering) > 0:
                        dis = abs(take_closest([it[0] for it in pair_after_filtering], p[0])-p[0])+abs(take_closest([it[1] for it in pair_after_filtering], p[1])-p[1])
                        if dis > self.MAX_TS_DIST and p[0]<j[0] and p[1]>j[-1]:
                            pair_after_filtering.append(p)
                    elif p[0]<j[0] and p[1]>j[-1]:
                        pair_after_filtering.append(p)
                    if len(pair_after_filtering) >= self.Max_site_per_splice: # only allow up to `Max_site_per_splice` combinations
                        break
            sup_cnt_total = 0
            for p in pair_after_filtering:
                sup_cnt_total += sum(it[1] for it in tmp_pair if abs(it[0][0]-p[0])+abs(it[0][1]-p[1])< self.MAX_TS_DIST)
            if float(sup_cnt_total)/len(self.lr_pair[j]) <= self.Min_sup_pct:
                # not enough support counts, often happens in the case when the TSS/TES is strongly degraded.
                #print "not enough support counts:",sup_cnt_total,len(self.lr_pair[j])
                if len(pair_enrich)==0:
                    for p, _ in tmp_pair:
                        pair_enrich.append((p, sum(it[1] for it in tmp_pair if abs(it[0][0]-p[0])+abs(it[0][1]-p[1])< self.MAX_TS_DIST)))
                    pair_enrich.sort(key=lambda x: x[1], reverse=True)
                tmp_ex = list(j)
                tmp_ex.insert(0,int(pair_enrich[0][0][0]))
                tmp_ex.append(int(pair_enrich[0][0][1]))
                self.raw_isoforms[tuple(tmp_ex)] = len(self.lr_pair[j])
                self.strand_cnt[tuple(tmp_ex)] = self.strand_cnt[j]
            else:
                for p in pair_after_filtering:  # add filtered TSS/TES to `raw_isoforms`
                    tmp_ex = list(j)
                    tmp_ex.insert(0,int(p[0]))
                    tmp_ex.append(int(p[1]))
                    self.raw_isoforms[tuple(tmp_ex)] = sum(it[1] for it in tmp_pair if abs(it[0][0]-p[0])+abs(it[0][1]-p[1])< self.MAX_TS_DIST)
                    self.strand_cnt[tuple(tmp_ex)] = self.strand_cnt[j]

    def match_known_annotation(self, transcript_to_junctions, transcript_dict, gene_dict, one_block, fa_dict):
        """
        annotate similar transcript and fix splicing position

        should run after `update_all_splice` amd `filter_TSS_TES`
        """
        Iso = namedtuple('Iso', ["support_cnt", "transcript_id","gene_id"])
        splice_site = get_splice_site(transcript_to_junctions, one_block.transcript_list)
        junc_list = [transcript_to_junctions[it]["junctions"] for it in one_block.transcript_list]
        junc_dict = dict((transcript_to_junctions[tr]["junctions"], tr) for tr in one_block.transcript_list)
        exons_list = [list(it) for it in junc_list]
        [exons_list[ith].append(transcript_to_junctions[tr]["right"]) for ith, tr in enumerate(one_block.transcript_list)]
        [exons_list[ith].insert(0,transcript_to_junctions[tr]["left"]) for ith, tr in enumerate(one_block.transcript_list)]
        exons_dict = dict((tuple(exons_list[ith]), tr) for ith, tr in enumerate(one_block.transcript_list))
        TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, one_block.transcript_list)
        if len(splice_site) == 0:
            return 0
        for one_exon in self.single_block_dict:
            for tr in one_block.transcript_list:
                if self.strand_specific==0:
                    tmp_std = self.strand_cnt[one_exon]
                else:
                    tmp_std = 1 if transcript_dict[tr].strand=="+" else -1
                if len(transcript_to_junctions[tr]["junctions"])==0 and (tmp_std == self.strand_cnt[one_exon]):
                    if abs(one_exon[0]-transcript_to_junctions[tr]["left"])<self.MAX_TS_DIST and abs(one_exon[1]-transcript_to_junctions[tr]["right"])<self.MAX_TS_DIST:
                        known_exons = (transcript_to_junctions[tr]["left"], transcript_to_junctions[tr]["right"])
                        if known_exons in self.known_isoforms:
                            self.known_isoforms[known_exons] = Iso(self.known_isoforms[known_exons].support_cnt+len(self.single_blocks[self.single_block_dict[one_exon]]),tr,transcript_dict[tr].parent_id)
                        else:
                            self.known_isoforms[known_exons] = Iso(len(self.single_blocks[self.single_block_dict[one_exon]]),tr,transcript_dict[tr].parent_id)
                            if self.strand_specific==0:  # if not strand specific protocol, use annotation
                                self.strand_cnt[tuple(known_exons)] = 1 if transcript_dict[tr].strand=="+" else -1
                            else:
                                self.strand_cnt[tuple(known_exons)] = self.strand_cnt[one_exon]
        for raw_iso in self.raw_isoforms:
            found = False
            if len(raw_iso)>len(set(raw_iso)):
                #print "REPEAT splice site, skip:", self.raw_isoforms[raw_iso], raw_iso
                continue
            for tr in one_block.transcript_list:
                if self.strand_specific==0:
                    tmp_std = self.strand_cnt[raw_iso]
                else:
                    tmp_std = 1 if transcript_dict[tr].strand=="+" else -1
                # same number of exons, same strand
                if (len(raw_iso)-2 == len(transcript_to_junctions[tr]["junctions"])) and (tmp_std == self.strand_cnt[raw_iso]):
                    if all(abs(it-ij)<self.MAX_DIST for it, ij in zip(raw_iso[1:-1], transcript_to_junctions[tr]["junctions"])):
                        if (abs(raw_iso[0]-transcript_to_junctions[tr]["left"]) < self.MAX_TS_DIST) and (abs(raw_iso[-1]-transcript_to_junctions[tr]["right"]) < self.MAX_TS_DIST):
                            known_exons = list(transcript_to_junctions[tr]["junctions"])
                            known_exons.insert(0,transcript_to_junctions[tr]["left"])
                            known_exons.append(transcript_to_junctions[tr]["right"])
                            if tuple(known_exons) in self.known_isoforms:
                                self.known_isoforms[tuple(known_exons)] = Iso(max(self.known_isoforms[tuple(known_exons)].support_cnt, self.raw_isoforms[raw_iso]),tr,transcript_dict[tr].parent_id)
                            else:
                                self.known_isoforms[tuple(known_exons)] = Iso(self.raw_isoforms[raw_iso],tr,transcript_dict[tr].parent_id)
                                if self.strand_specific==0:  # if not strand specific protocol, use annotation
                                    self.strand_cnt[tuple(known_exons)] = 1 if transcript_dict[tr].strand=="+" else -1
                                else:
                                    self.strand_cnt[tuple(known_exons)] = self.strand_cnt[raw_iso]
                            found = True
                            break
            if not found:
                new_exons = find_best_splice_chain(raw_iso, junc_list, self.MAX_SPLICE_MATCH_DIST)
                if not all(new_exons[ix]<new_exons[ix+1] for ix in range(len(new_exons)-1)):
                    new_exons = list(raw_iso)
                for ix, a_site in enumerate(new_exons):
                    if ix==0:  # left
                        clo = take_closest(TSS_TES_site["left"],a_site)
                        if abs(clo-a_site)<self.MAX_TS_DIST and clo < raw_iso[ix+1]:
                            new_exons[ix] = clo
                        else:
                            new_exons[ix] = a_site
                    elif 0<ix<len(raw_iso)-1:
                        clo = take_closest(splice_site,a_site)
                        if abs(clo-a_site)<self.MAX_SPLICE_MATCH_DIST and clo > new_exons[-1] and clo < raw_iso[ix+1]:
                            new_exons[ix] = clo
                        else:
                            new_exons[ix] = a_site
                    else:  # right
                        clo = take_closest(TSS_TES_site["right"],a_site)
                        if abs(clo-a_site)<self.MAX_TS_DIST and clo > new_exons[-1]:
                            new_exons[ix] = clo
                        else:
                            new_exons[ix] = a_site
                if all(new_exons[ix]<new_exons[ix+1] for ix in range(len(new_exons)-1)):
                    #if tuple(new_exons[1:-1]) in junc_dict and abs(new_exons[0]-transcript_to_junctions[junc_dict[tuple(new_exons[1:-1])]]["left"])<self.MAX_TS_DIST and abs(new_exons[-1]-transcript_to_junctions[junc_dict[tuple(new_exons[1:-1])]]["right"])<self.MAX_TS_DIST:
                    #    known_exons = list(transcript_to_junctions[ junc_dict[tuple(new_exons[1:-1])] ]["junctions"])
                    #    known_exons.insert(0,transcript_to_junctions[ junc_dict[tuple(new_exons[1:-1])] ]["left"])
                    #    known_exons.append(transcript_to_junctions[ junc_dict[tuple(new_exons[1:-1])] ]["right"])
                    #    if tuple(known_exons) not in self.known_isoforms:
                    #        self.known_isoforms[tuple(known_exons)] = Iso(self.raw_isoforms[raw_iso],junc_dict[tuple(new_exons[1:-1])],transcript_dict[junc_dict[tuple(new_exons[1:-1])]].parent_id)
                    #        self.strand_cnt[tuple(known_exons)] = 1 if transcript_dict[junc_dict[tuple(new_exons[1:-1])]].strand=="+" else -1
                    if tuple(new_exons) not in self.new_isoforms:
                        self.new_isoforms[tuple(new_exons)] = Iso(self.raw_isoforms[raw_iso],"","")
                        self.strand_cnt[tuple(new_exons)] = self.strand_cnt[raw_iso]
                    else:
                        self.new_isoforms[tuple(new_exons)] = Iso(self.new_isoforms[tuple(new_exons)].support_cnt+self.raw_isoforms[raw_iso],"","")
                else:
                    print "exon chain not in ascending order:", raw_iso,new_exons
        
        # remove incomplete transcript (due to 3' bias)
        if self.remove_incomp_reads>0:
            del_key = []
            for ni in self.new_isoforms:  # 1. match to known isoform detected and remove new isoform if within known isoform
                if ni in self.known_isoforms:
                    del_key.append(ni)
                    continue
                for kn in self.known_isoforms:
                    if len(ni) < len(kn) and if_exon_contains(kn, ni[2:], 1):
                        s_l = fa_dict[self.ch][(ni[0]-15):ni[0]]
                        s_r = fa_dict[self.ch][ni[0]:(ni[0]+15)]
                        if s_l.count('T')>10 or s_l.count('A')>10 or s_r.count('T')>10 or s_r.count('A')>10:
                            del_key.append(ni)
                            break
                    elif len(ni) < len(kn) and if_exon_contains(kn, ni[:-2], 1):
                        s_l = fa_dict[self.ch][(ni[-1]-15):ni[-1]]
                        s_r = fa_dict[self.ch][ni[-1]:(ni[-1]+15)]
                        if s_l.count('T')>10 or s_l.count('A')>10 or s_r.count('T')>10 or s_r.count('A')>10:
                            del_key.append(ni)
                            break
                    if len(ni) < len(kn) and if_exon_contains(kn, ni, self.MAX_TS_DIST):
                        sim_pct_sq = get_exon_sim_pct(ni, kn)**2
                        if  (self.new_isoforms[ni].support_cnt < (1+sim_pct_sq*self.remove_incomp_reads)*self.known_isoforms[kn].support_cnt):
                            del_key.append(ni)
                            break
                    if ni==kn:
                        del_key.append(ni)  # should not be the same
                        break
            for key in list(set(del_key)):
                del self.new_isoforms[key]
            del_key = []
            for ni in self.new_isoforms:  # 2. match to known isoform in annotation and remove new isoform if very similar to annotation
                for kn in exons_dict:
                    if kn not in self.known_isoforms:
                        if len(ni) < len(kn) and if_exon_contains(kn, ni, self.MAX_TS_DIST):
                            sim_pct = get_exon_sim_pct(ni, kn)
                            if sim_pct>0.95:
                                del_key.append(ni)
                                self.known_isoforms[kn] = Iso(self.new_isoforms[ni].support_cnt,exons_dict[kn],transcript_dict[exons_dict[kn]].parent_id)
                                self.strand_cnt[kn] = 1 if transcript_dict[exons_dict[kn]].strand=="+" else -1
                        elif ni==kn:
                            del_key.append(ni)  # should not be the same
                            self.known_isoforms[kn] = Iso(self.new_isoforms[ni].support_cnt,exons_dict[kn],transcript_dict[exons_dict[kn]].parent_id)
                            self.strand_cnt[kn] = 1 if transcript_dict[exons_dict[kn]].strand=="+" else -1
            for key in list(set(del_key)):
                del self.new_isoforms[key]
            if len(self.new_isoforms) > 1:  # 3. match among new isoform
                del_key = []
                isos = self.new_isoforms.keys()
                for i in range(len(isos)-1):
                    for j in range(i+1,len(isos)):
                        if len(isos[i]) < len(isos[j]) and if_exon_contains(isos[j],isos[i],self.MAX_TS_DIST):
                            sim_pct_sq = get_exon_sim_pct(isos[j], isos[i])**2
                            if (self.new_isoforms[isos[i]].support_cnt < sim_pct_sq*self.remove_incomp_reads*self.new_isoforms[isos[j]].support_cnt):
                                del_key.append(isos[i])
                        elif len(isos[i]) > len(isos[j]) and if_exon_contains(isos[i],isos[j],self.MAX_TS_DIST):
                            sim_pct_sq = get_exon_sim_pct(isos[j], isos[i])**2
                            if (self.new_isoforms[isos[j]].support_cnt < sim_pct_sq*self.remove_incomp_reads*self.new_isoforms[isos[i]].support_cnt):
                                del_key.append(isos[j])
                for key in list(set(del_key)):
                    del self.new_isoforms[key]


        # remove reads start or end with polyA (internal priming)
        #del_key = []
        #for ni in self.new_isoforms:
        #    s_l = fa_dict[self.ch][(ni[-1]-15):ni[0]]
        #    s_r = fa_dict[self.ch][ni[-1]:(ni[-1]+15)]
        #    if s_l.count('T')>10 or s_l.count('A')>10 or s_r.count('T')>10 or s_r.count('A')>10:
        #        del_key.append(ni)
        #for key in del_key:
        #    del self.new_isoforms[key]


        # match to gene
        del_key = []
        update_iso_dict = {}
        if len(self.new_isoforms)>0:
            for i in self.new_isoforms:
                iso_len = sum(it[1]-it[0] for it in pairwise(i))
                tmp = [(exon_overlap(i, gene_dict[ge]),ge) for ge in one_block.gene_to_tr]
                if self.strand_specific != 0:  # if  strand specific protocol, use read
                    stnd = "+" if self.strand_cnt[i]==1 else "-"
                    tmp = [it for it in tmp if transcript_dict[one_block.gene_to_tr[it[1]][0]].strand == stnd]
                    if len(tmp) == 0:
                        #print "no match", stnd, one_block.gene_to_tr.keys()
                        continue
                tmp.sort(key=lambda x:x[0],reverse = True)
                if tmp[0][0]>0:
                    if i[0]>= gene_dict[tmp[0][1]][0] and i[-1]<= gene_dict[tmp[0][1]][-1]:
                        self.new_isoforms[i] = Iso(self.new_isoforms[i].support_cnt, "", tmp[0][1])
                    elif exon_overlap(i[2:4], gene_dict[tmp[0][1]])>0 and exon_overlap(i[-4:-2], gene_dict[tmp[0][1]])>0:  # the second and the second last within the gene
                        if i[1]-i[0]<=self.min_fl_exon_len and exon_overlap(i[0:2], gene_dict[tmp[0][1]])==0:
                            ba = fa_dict[self.ch][i[0]:i[1]]
                            if ba.count('T') > 0.7*len(ba) or ba.count('A') > 0.7*len(ba):  # leftover polyA tail
                                del_key.append(i)
                                update_iso_dict[i[2:]] = Iso(self.new_isoforms[i].support_cnt, "", tmp[0][1])
                        elif i[-1]-i[-2]<=self.min_fl_exon_len and exon_overlap(i[-2:], gene_dict[tmp[0][1]])==0:
                            ba = fa_dict[self.ch][i[-2]:i[-1]]
                            if ba.count('T') > 0.7*len(ba) or ba.count('A') > 0.7*len(ba):  # leftover polyA tail
                                del_key.append(i)
                                update_iso_dict[i[:-2]] = Iso(self.new_isoforms[i].support_cnt, "", tmp[0][1])
                        else:
                            self.new_isoforms[i] = Iso(self.new_isoforms[i].support_cnt, "", tmp[0][1])
                    else:
                        if tmp[0][0] > 0.8*iso_len:
                            if (i[1]-i[0]<self.min_fl_exon_len) or (i[-1]-i[-2]<self.min_fl_exon_len):
                                continue  # alignment artifact
                            else:
                                self.new_isoforms[i] = Iso(self.new_isoforms[i].support_cnt, "", tmp[0][1])  # might be real eRNA
                        elif sum(1 if it in splice_site else 0 for it in i)>4:
                            self.new_isoforms[i] = Iso(self.new_isoforms[i].support_cnt, "", tmp[0][1])  # more than 4 splice site match
                        else:
                            pass
                            #print self.ch, i
                    # looking at most overlaped gene, if no overlap to any gene, keep empty
                    # TODO: better at fusion genes
            if len(del_key)>0:
                for key in del_key:
                    del self.new_isoforms[key]
                #for u_i in update_iso_dict:
                #    self.new_isoforms[u_i] = update_iso_dict[u_i]
                #    self.strand_cnt[u_i] = self.strand_cnt[del_key[0]]
            for i in self.new_isoforms:
                if self.new_isoforms[i].gene_id == "":
                    # print "no matching genes:"
                    # print i, self.new_isoforms[i].support_cnt
                    continue
                self.ge_dict.setdefault(self.new_isoforms[i].gene_id, []).append(i)
                if self.strand_specific == 0:
                    if self.strand_cnt[i] == 0:
                        self.strand_cnt[i] = 1 if transcript_dict[one_block.gene_to_tr[self.new_isoforms[i].gene_id][0]].strand=="+" else -1
        for i in self.known_isoforms:
            self.ge_dict.setdefault(self.known_isoforms[i].gene_id, []).append(i)
                        
    def isoform_to_gff3(self, isoform_pct=-1):
        gff3_fmt = "{_ch}\t{_sr}\t{_ty}\t{_st}\t{_en}\t{_sc}\t{_stnd}\t{_ph}\t{_attr}"
        gff_rec = []
        transcript_id_dict={}
        if len(self.new_isoforms)+len(self.known_isoforms) == 0:
            return ""
        for g in self.ge_dict:
            gff_tmp = []
            total_cnt = sum(self.new_isoforms[e].support_cnt for e in self.ge_dict[g] if e in self.new_isoforms)+\
            sum(self.known_isoforms[e].support_cnt for e in self.ge_dict[g] if e in self.known_isoforms)
            gff_tmp.append(gff3_fmt.format(_ch=self.ch,_sr=".",_ty="gene",
            _st=min(ix[0] for ix in self.ge_dict[g] if ix in self.known_isoforms or self.new_isoforms[ix].support_cnt>isoform_pct*total_cnt)+1,
            _en=max(ix[-1] for ix in self.ge_dict[g] if ix in self.known_isoforms or self.new_isoforms[ix].support_cnt>isoform_pct*total_cnt),
            _sc=".",_stnd="+" if self.strand_cnt[self.ge_dict[g][0]]==1 else "-", _ph=".",
            _attr="ID=gene:{};gene_id={};support_count={}".format(g, g, total_cnt)))
            for exons in self.ge_dict[g]:
                if exons in self.new_isoforms and exons in self.known_isoforms:
                    print "BOTH in new and known:",self.new_isoforms[exons],self.known_isoforms[exons]
                if exons in self.new_isoforms:
                    source = "new"
                    sup_cnt = self.new_isoforms[exons].support_cnt
                    if 0<isoform_pct<1 and sup_cnt<(isoform_pct*total_cnt):
                        continue
                    tp_id = "{}_{}_{}".format(g,exons[0]+1,exons[-1])
                    if tp_id not in transcript_id_dict:
                        transcript_id_dict[tp_id] = 1
                        tp_id = "{}_1".format(tp_id)
                    else:
                        transcript_id_dict[tp_id] += 1
                        tp_id = "{}_{}".format(tp_id, transcript_id_dict[tp_id])
                else:
                    source = "known"
                    sup_cnt = self.known_isoforms[exons].support_cnt
                    tp_id = self.known_isoforms[exons].transcript_id
                if sup_cnt<self.Min_sup_cnt:
                    continue
                exon_idx = 1
                gff_tmp.append(gff3_fmt.format(_ch=self.ch,_sr=source,_ty="transcript",_st=exons[0]+1,_en=exons[-1],
                            _sc=".",_stnd="+" if self.strand_cnt[exons]==1 else "-", _ph=".",
                            _attr="ID=transcript:{};transcript_id={};Parent=gene:{};support_count={};source={}".format(tp_id, tp_id, g, sup_cnt, source )))
                for ex in range(0,len(exons),2):
                    gff_tmp.append(gff3_fmt.format(_ch=self.ch,_sr=source,_ty="exon",_st=exons[ex]+1,_en=exons[ex+1], # `+1` because gff is 1-based
                        _sc=".",_stnd="+" if self.strand_cnt[exons]==1 else "-", _ph=".",
                        _attr="exon_id=exon:{}_{};Parent=transcript:{};rank={}".format(exons[ex]+1, exons[ex+1],tp_id, exon_idx )))
                    exon_idx += 1
            if len(gff_tmp) >2:
                gff_rec.extend(gff_tmp)
        if len(gff_rec)>0:
            return "\n".join(gff_rec)+"\n"
        else:
            return ""

    def raw_splice_to_gff3(self):
        """
        write all raw splice form, only for debug purpose
        """
        gff3_fmt = "{_ch}\t{_sr}\t{_ty}\t{_st}\t{_en}\t{_sc}\t{_stnd}\t{_ph}\t{_attr}"
        gff_tmp = []
        transcript_id_dict={}
        exon_id_dict={}
        for iso in self.raw_isoforms:
            exon_idx = 1
            tp_id = "{}_{}".format(iso[0]+1,iso[-1])
            if tp_id not in transcript_id_dict:
                transcript_id_dict[tp_id] = 1
                tp_id = "{}_1".format(tp_id)
            else:
                transcript_id_dict[tp_id] += 1
                tp_id = "{}_{}".format(tp_id, transcript_id_dict[tp_id])
            # transcript
            gff_tmp.append(gff3_fmt.format(_ch=self.ch,_sr="raw",_ty="transcript",_st=iso[0]+1,_en=iso[-1],
                _sc=".",_stnd="+" if self.strand_cnt[iso]==1 else "-", _ph=".",
                _attr="ID=transcript:{};support_count={}".format(tp_id, self.raw_isoforms[iso] )))
            for ex in range(0,len(iso),2):
                gff_tmp.append(gff3_fmt.format(_ch=self.ch,_sr="raw",_ty="exon",_st=iso[ex]+1,_en=iso[ex+1], # `+1` because gff is 1-based
                    _sc=".",_stnd="+" if self.strand_cnt[iso]==1 else "-", _ph=".",
                    _attr="exon_id=exon:{}_{};Parent=transcript:{};rank={}".format(iso[ex]+1, iso[ex+1],tp_id, exon_idx )))
                exon_idx += 1
        if len(gff_tmp)>0:
            return "\n".join(gff_tmp)+"\n"
        else:
            return ""





def group_bam2isoform(bam_in, out_gff3, out_stat, summary_csv, chr_to_blocks, gene_dict, transcript_to_junctions, transcript_dict, fa_f, config, downsample_ratio, raw_gff3=None):
    if "random_seed" in config.keys():
        random.seed(config["random_seed"])
    else:
        random.seed(666666)
    bamfile = pysam.AlignmentFile(bam_in, "rb")
    #csv_out = open(summary_csv,"w")
    iso_annotated = open(out_gff3,"w")
    iso_annotated.write("##gff-version 3\n")
    if raw_gff3 is not None:
        splice_raw = open(raw_gff3,"w")
        splice_raw.write("##gff-version 3\n")
    tss_tes_stat = open(out_stat,"w")
    isoform_dict = {}
    fa_dict = {}
    for c in get_fa(fa_f):
        fa_dict[c[0]] = c[1]
    for ch in chr_to_blocks:
        # print ch
        #if ch != "5":
        #    continue
        print ch
        for ith, bl in enumerate(chr_to_blocks[ch]):
            it_region = bamfile.fetch(ch, bl.s, bl.e)
            TSS_TES_site = get_TSS_TES_site(transcript_to_junctions, bl.transcript_list)
            tmp_isoform = Isoforms(ch, config)
            for rec in it_region:
                if 0<downsample_ratio<1 and random.uniform(0, 1)>downsample_ratio:
                    continue   # downsample analysis
                if rec.is_secondary:
                    continue
                rec.cigar = smooth_cigar(rec.cigar, thr=20)
                blocks = rec.get_blocks()
                junctions = blocks_to_junctions(blocks)
                tmp_isoform.add_isoform(junctions,rec.is_reverse)
            if len(tmp_isoform)>0:
                tmp_isoform.update_all_splice()
                tmp_isoform.filter_TSS_TES(tss_tes_stat,known_site=TSS_TES_site,fdr_cutoff=0.1)
                #tmp_isoform.site_stat(tss_tes_stat)
                tmp_isoform.match_known_annotation(transcript_to_junctions, transcript_dict, gene_dict, bl, fa_dict)
                isoform_dict[(ch, bl.s, bl.e)] =tmp_isoform
                if raw_gff3 is not None:
                    splice_raw.write(tmp_isoform.raw_splice_to_gff3())
                iso_annotated.write(tmp_isoform.isoform_to_gff3(isoform_pct=config["Min_cnt_pct"]))
    #with open(iso_exact,"w") as out_f:
    #    out_f.write("##gff-version 3\n")
    tss_tes_stat.close()
    iso_annotated.close()
    bamfile.close()
    if raw_gff3 is not None:
        splice_raw.close()





if __name__ == '__main__':
    pass

