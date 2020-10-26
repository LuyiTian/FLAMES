# filter gff file after realign to transcript.
from parse_gene_anno import parse_gff_tree
from collections import Counter
import copy

def annotate_full_splice_match(transcript_to_exon,transcript_to_exon_ref,transcript_dict,transcript_dict_ref,anno_out,tr_cnt,min_sup_reads):
    splice_dict = {}
    splice_dict_ref = {}
    splice_gene_dict = {}
    tr_id_list = []
    FSM_id_list = []
    FSM_to_reference = []
    gene_id_list = []
    for tr in transcript_to_exon:
        tmp = []
        for ex in transcript_to_exon[tr]:
            tmp.append(ex[0])
            tmp.append(ex[1])
        if len(tmp)==2:
            continue
        tmp = tuple(tmp[1:-1])
        if (tr in tr_cnt) and tr_cnt[tr]>=min_sup_reads:
            splice_dict.setdefault(tmp,[]).append(tr)
            splice_gene_dict.setdefault(tmp,[]).append(transcript_dict[tr].parent_id)
    for tr in transcript_to_exon_ref:
        tmp = []
        for ex in transcript_to_exon_ref[tr]:
            tmp.append(ex[0])
            tmp.append(ex[1])
        if len(tmp)==2:
            continue
        tmp = tuple(tmp[1:-1])
        if (tr in tr_cnt) and (tr not in transcript_to_exon) and (tr_cnt[tr]>=min_sup_reads):
            splice_dict.setdefault(tmp,[]).append(tr)
            splice_gene_dict.setdefault(tmp,[]).append(transcript_dict_ref[tr].parent_id)
        splice_dict_ref.setdefault(tmp,[]).append(tr)
    for s in splice_gene_dict:
        if len(set(splice_gene_dict[s]))>1:
            print "shared splice chain:",set(splice_gene_dict[s])
            print s
    for s in splice_dict:
        if len(splice_dict[s])==1:
            tr_id_list.append(splice_dict[s][0])
            gene_id_list.append(splice_gene_dict[s][0])
            if s in splice_dict_ref:
                FSM_id_list.append(splice_dict_ref[s][0])
                FSM_to_reference.append(True)
            else:
                FSM_id_list.append(splice_dict[s][0])
                FSM_to_reference.append(False)
        else:
            ref_tr = [tr for tr in splice_dict[s] if tr in transcript_to_exon_ref]
            if len(ref_tr)>0:
                for tr in splice_dict[s]:
                    tr_id_list.append(tr)
                    gene_id_list.append(splice_gene_dict[s][0])
                    FSM_id_list.append(ref_tr[0])
                    FSM_to_reference.append(True)
            elif s in splice_dict_ref:
                for tr in splice_dict[s]:
                    tr_id_list.append(tr)
                    gene_id_list.append(splice_gene_dict[s][0])
                    FSM_id_list.append(splice_dict_ref[s][0])
                    FSM_to_reference.append(True)
            else:
                for tr in splice_dict[s]:
                    tr_id_list.append(tr)
                    gene_id_list.append(splice_gene_dict[s][0])
                    FSM_id_list.append(splice_dict[s][0])
                    FSM_to_reference.append(False)
    with open(anno_out,"w") as f:
        f.write("transcript_id,gene_id,FSM_match,FSM_match_to_ref,total_count\n")
        for i in range(len(tr_id_list)):
            f.write("{},{},{},{},{}\n".format(tr_id_list[i],gene_id_list[i],FSM_id_list[i],FSM_to_reference[i],tr_cnt[tr_id_list[i]]))


def annotate_filter_gff(isoform_gff,ref_gff,isoform_out,anno_out,tr_cnt,min_sup_reads,verbose=True):
    """
    combine FLAMES ouput with reference and filter out transcript by
    realignment result
    """
    gff3_fmt = "{_ch}\t{_sr}\t{_ty}\t{_st}\t{_en}\t{_sc}\t{_stnd}\t{_ph}\t{_attr}"

    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(isoform_gff)
    chr_to_gene_ref, transcript_dict_ref, gene_to_transcript_ref, transcript_to_exon_ref = parse_gff_tree(ref_gff)
    prt = "Filtering and combining isoforms from realigned bam file:\n\tBefore filtering: {} isoforms in count matrix. {} isoforms in reference annotation. {} isoforms in FLAMES raw output.".format(
        len(tr_cnt),
        len(transcript_to_exon_ref),
        len(transcript_to_exon)
    )
    if verbose:
        print(prt)
    gff_rec = []
    iso_rm=0
    iso_kp=0
    for ch in chr_to_gene:
        new_ge_list = copy.deepcopy(chr_to_gene[ch])
        new_ge_list.extend(chr_to_gene_ref[ch])
        new_ge_list = list(set(new_ge_list))
        for ge in new_ge_list:
            gff_tmp = []
            total_cnt = 0
            mi = 99999999999
            ma = -1
            if ge in gene_to_transcript:
                for tr in gene_to_transcript[ge]:
                    if (tr in tr_cnt) and (tr_cnt[tr]>=min_sup_reads):
                        gff_tmp.append(gff3_fmt.format(_ch=ch,_sr="FLAMES",_ty="transcript",
                            _st=transcript_to_exon[tr][0][0]+1,
                            _en=transcript_to_exon[tr][-1][1],
                            _sc=".",_stnd=transcript_dict[tr].strand, _ph=".",
                            _attr="ID=transcript:{};transcript_id={};Parent=gene:{};support_count={}".format(tr,tr, ge, tr_cnt[tr])))
                        total_cnt += tr_cnt[tr]
                        mi = min(mi,transcript_to_exon[tr][0][0])
                        ma = max(ma,transcript_to_exon[tr][-1][1])
                        exon_idx = 1
                        for ex in transcript_to_exon[tr]:
                            gff_tmp.append(gff3_fmt.format(_ch=ch,_sr="FLAMES",_ty="exon",_st=ex[0]+1,_en=ex[1], # `+1` because gff is 1-based
                            _sc=".",_stnd=transcript_dict[tr].strand, _ph=".",
                            _attr="exon_id=exon:{}_{};Parent=transcript:{};rank={}".format(ex[0]+1, ex[1],tr, exon_idx )))
                            exon_idx += 1
                        iso_kp += 1
                    else:
                        iso_rm += 1
            if ge in gene_to_transcript_ref:
                for tr in gene_to_transcript_ref[ge]:
                    if (tr not in transcript_dict) and (tr in tr_cnt) and (tr_cnt[tr]>=min_sup_reads):  # not in FLAMES output but in tr count
                        gff_tmp.append(gff3_fmt.format(_ch=ch,_sr="reference",_ty="transcript",
                            _st=transcript_to_exon_ref[tr][0][0]+1,
                            _en=transcript_to_exon_ref[tr][-1][1],
                            _sc=".",_stnd=transcript_dict_ref[tr].strand, _ph=".",
                            _attr="ID=transcript:{};transcript_id={};Parent=gene:{};support_count={}".format(tr,tr, ge, tr_cnt[tr])))
                        total_cnt += tr_cnt[tr]
                        mi = min(mi,transcript_to_exon_ref[tr][0][0])
                        ma = max(ma,transcript_to_exon_ref[tr][-1][1])
                        exon_idx = 1
                        for ex in transcript_to_exon_ref[tr]:
                            gff_tmp.append(gff3_fmt.format(_ch=ch,_sr="reference",_ty="exon",_st=ex[0]+1,_en=ex[1], # `+1` because gff is 1-based
                            _sc=".",_stnd=transcript_dict_ref[tr].strand, _ph=".",
                            _attr="exon_id=exon:{}_{};Parent=transcript:{};rank={}".format(ex[0]+1, ex[1],tr, exon_idx )))
                            exon_idx += 1
                        iso_kp += 1
            if len(gff_tmp)>0 and ge in gene_to_transcript:
                gff_tmp.insert(0,gff3_fmt.format(_ch=ch,_sr="FLAMES",_ty="gene",
                    _st=mi+1,
                    _en=ma,
                    _sc=".",_stnd=transcript_dict[gene_to_transcript[ge][0]].strand, _ph=".",
                    _attr="ID=gene:{};gene_id={};support_count={}".format(ge, ge, total_cnt)))
                gff_rec.extend(gff_tmp)
            elif len(gff_tmp)>0:
                gff_tmp.insert(0,gff3_fmt.format(_ch=ch,_sr="FLAMES",_ty="gene",
                    _st=mi+1,
                    _en=ma,
                    _sc=".",_stnd=transcript_dict_ref[gene_to_transcript_ref[ge][0]].strand, _ph=".",
                    _attr="ID=gene:{};gene_id={};support_count={}".format(ge, ge, total_cnt)))
                gff_rec.extend(gff_tmp)
    iso_annotated = open(isoform_out,"w")
    iso_annotated.write("##gff-version 3\n")
    iso_annotated.write("\n".join(gff_rec))
    iso_annotated.close()
    prt = "\tAfter filtering: kept {} isoforms. removed {} isoforms.".format(
        iso_kp,
        iso_rm
    )
    if verbose:
        print(prt)
    annotate_full_splice_match(transcript_to_exon,transcript_to_exon_ref,transcript_dict,transcript_dict_ref,anno_out,tr_cnt,min_sup_reads)