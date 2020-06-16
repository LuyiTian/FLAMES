# create merged fastq file for bulk long read RNAseq
from collections import Counter
import gzip
import os
import random
import string

class Fastq(object):
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]

    def __repr__(self):
        return "Fastq({name})".format(name=self.name)

    def __str__(self):
        return "@{name}\n{seq}\n+\n{qual}\n".format(
            name=self.name, seq=self.seq, qual=self.qual)
    def to_fa(self):
        return ">{name}\n{seq}\n".format(name=self.name, seq=self.seq)


def readfq(fq):
    record = []
    openFunc = gzip.open if fq.endswith(".gz") else open
    for line in openFunc(fq):
        record.append(line.strip())
        if len(record) == 4:
            yield Fastq(record)
            record = []


def merge_bulk_fq(fq_dir, anno_csv, out_fq):
    fq_names = [it for it in os.listdir(fq_dir) if ("fq" in it) or ("fastq" in it)]
    fq_dict = {}
    for f in fq_names:
        fq_dict[f] = os.path.join(fq_dir,f)
    merged_fq = gzip.open(out_fq, 'wb')
    pseudo_bc_dict = {}
    fq_cnt = Counter()
    random.seed(2333666)
    for a_fq in fq_dict:
        while True:
            pseudo_bc = ''.join(random.choice(string.ascii_uppercase) for _ in range(16))
            if pseudo_bc not in pseudo_bc_dict:
                print a_fq, pseudo_bc
                break
        pseudo_bc_dict[pseudo_bc] = a_fq
        for rec in readfq(fq_dict[a_fq]):
            rec.name = "{}_NNN#{}".format(pseudo_bc, rec.name)
            merged_fq.write(rec.__str__())
            fq_cnt[a_fq] += 1
    merged_fq.close()
    with open(anno_csv,"w") as f:
        f.write("file_name,pseudo_barcode\n")
        for pseudo_bc in pseudo_bc_dict:
            f.write("{},{}\n".format(pseudo_bc_dict[pseudo_bc],pseudo_bc))
    for i in fq_cnt:
        print i, fq_cnt[i]
