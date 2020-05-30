# gff3 to fasta

from parse_gene_anno import parse_gff_tree
import subprocess
CP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N","a":"t","t":"a","c":"g","g":"c"}


def r_c(seq):
    new_seq = []
    for c in seq[::-1]:
        new_seq.append(CP[c])
    return "".join(new_seq)


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
            seq.append(line.strip())
    yield ch, "".join(seq)


def write_fa(fn, na, seq, warp_len=50):
    fn.write(">{}\n".format(na))
    ix = 0
    while ix < len(seq):
        if ix+warp_len > len(seq):
            fn.write("{}\n".format(seq[ix:(len(seq))]))
            ix = ix+warp_len
        else:
            fn.write("{}\n".format(seq[ix:(ix+warp_len)]))
            ix = ix+warp_len


def get_transcript_seq(fa_file, fa_out_f, chr_to_gene, transcript_dict,
                       gene_to_transcript, transcript_to_exon, ref_dict = None):
    global_isoform_dict = {}
    global_seq_dict = {}
    fa_dict = {}
    fa_out = open(fa_out_f, "w")
    for ch, seq in get_fa(fa_file):
        if ch not in chr_to_gene:
            continue
        else:
            pass
            # print("start to process chromosome", ch)
        for gene in chr_to_gene[ch]:
            for tr in gene_to_transcript[gene]:
                iso_l = []
                for e in transcript_to_exon[tr]:
                    assert (e[0]<e[1]),"exon end should be greater than exon start position."
                    iso_l.append(e[0])
                    iso_l.append(e[1])
                if tuple(iso_l) in global_isoform_dict:
                    print "duplicate transcript annotation:", global_isoform_dict[tuple(iso_l)], tr
                else:
                    global_isoform_dict[tuple(iso_l)] = tr
                    tr_seq = []
                    for e in transcript_to_exon[tr]:
                        tr_seq.append(seq[e[0]:e[1]])
                    tr_seq = "".join(tr_seq)
                    if transcript_dict[tr].strand != "+":
                        tr_seq = r_c(tr_seq)
                    fa_dict[tr] = tr_seq
                    if tr_seq in global_seq_dict:
                        print "duplicate transcript sequence:", global_seq_dict[tr_seq], tr
                    else:
                        global_seq_dict[tr_seq] = tr
                        #write_fa(fa_out, tr, tr_seq)
            if ref_dict is not None:
                if gene in ref_dict["chr_to_gene"][ch]:
                    for tr in ref_dict["gene_to_transcript"][gene]:
                        if tr in transcript_to_exon:
                            continue  # already in new annotation
                        iso_l = []
                        for e in ref_dict["transcript_to_exon"][tr]:
                            assert (e[0]<e[1]),"exon end should be greater than exon start position."
                            iso_l.append(e[0])
                            iso_l.append(e[1])
                        if tuple(iso_l) in global_isoform_dict:
                            if global_isoform_dict[tuple(iso_l)] not in ref_dict["transcript_dict"]:
                                print "transcript with same coordination", global_isoform_dict[tuple(iso_l)],tr
                                print iso_l
                                global_seq_dict[fa_dict[global_isoform_dict[tuple(iso_l)]]] = tr
                        else:
                            global_isoform_dict[tuple(iso_l)] = tr
                            tr_seq = []
                            for e in ref_dict["transcript_to_exon"][tr]:
                                tr_seq.append(seq[e[0]:e[1]])
                            tr_seq = "".join(tr_seq)
                            if ref_dict["transcript_dict"][tr].strand != "+":
                                tr_seq = r_c(tr_seq)
                            if tr_seq in global_seq_dict:
                                print "duplicate transcript sequence:", global_seq_dict[tr_seq], tr
                                global_seq_dict[tr_seq] = tr
                            else:
                                global_seq_dict[tr_seq] = tr
                                #write_fa(fa_out, tr, tr_seq)
    for tr_seq in global_seq_dict:
        write_fa(fa_out, global_seq_dict[tr_seq], tr_seq)
    fa_out.close()
    print subprocess.check_output(["samtools faidx {}".format(fa_out_f)], shell=True, stderr=subprocess.STDOUT)





if __name__ == '__main__':
    gff_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoforms/isoform_annotated.sample.nofilter.gff3"
    fa_file = "/stornext/General/data/user_managed/grpu_mritchie_1/LuyiTian/Index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    fa_out_f = "/stornext/General/data/user_managed/grpu_mritchie_1/SCmixology/PromethION/isoforms/human_GRCh38_transcript.sample.fa"
    chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon = parse_gff_tree(gff_f)
    get_transcript_seq(fa_file, fa_out_f, chr_to_gene, transcript_dict, gene_to_transcript, transcript_to_exon)
