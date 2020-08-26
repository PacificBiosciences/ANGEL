#!/usr/bin/env python
import os, sys
from collections import defaultdict
from Bio import SeqIO
__author__ = 'etseng@pacb.com'

def compare_files(input_prefix1, input_prefix2, output_prefix):
    cds1 = input_prefix1 + ".cds"
    utr1 = input_prefix1 + ".utr"
    pep1 = input_prefix1 + ".pep"
    cds2 = input_prefix2 + ".cds"
    utr2 = input_prefix2 + ".utr"
    pep2 = input_prefix2 + ".pep"

    if not os.path.exists(cds1) or not os.path.exists(utr1) or not os.path.exists(pep1):
        print("Abort! One or more of following files missing:", cds1, utr1, pep1, file=sys.stderr)
    if not os.path.exists(cds2) or not os.path.exists(utr2) or not os.path.exists(pep2):
        print("Abort! One or more of following files missing:", cds2, utr2, pep2, file=sys.stderr)

    d1 = defaultdict(lambda: 0)
    for r in SeqIO.parse(open(pep1), 'fasta'):
        _id = r.id.split('|')[0] # ex: PB.3.1
        d1[_id] = max(len(r.seq), d1[_id])

    d2 = defaultdict(lambda: 0)
    for r in SeqIO.parse(open(pep2), 'fasta'):
        _id = r.id.split('|')[0] # ex: PB.3.1
        d2[_id] = max(len(r.seq), d2[_id])

    f_out = open(output_prefix+'.compare.txt', 'w')
    f_out.write("#file1:{0}\n".format(input_prefix1))
    f_out.write("#file2:{}\n".format(input_prefix2))
    f_out.write("id\tlen1\tlen2\tpick\n")

    to_use = {} # seq id --> 1 if to use file1, 2 otherwise
    keys = list(set(d1.keys()).union(list(d2.keys())))
    keys.sort()
    for k in keys:
        if d1[k] >= d2[k]:
            to_use[k] = 1
        else:
            to_use[k] = 2
        f_out.write("{0}\t{1}\t{2}\t{3}\n".format(k, d1[k], d2[k], to_use[k]))
    f_out.close()


    f_cds = open(output_prefix+'.cds', 'w')
    f_pep = open(output_prefix+'.pep', 'w')
    f_utr = open(output_prefix+'.utr', 'w')


    for r in SeqIO.parse(open(cds1), 'fasta'):
        _id = r.id.split('|')[0]
        if to_use[_id] == 1:
            SeqIO.write(r, f_cds, 'fasta')
    for r in SeqIO.parse(open(pep1), 'fasta'):
        _id = r.id.split('|')[0]
        if to_use[_id] == 1:
            SeqIO.write(r, f_pep, 'fasta')
    for r in SeqIO.parse(open(utr1), 'fasta'):
        _id = r.id.split('|')[0]
        if to_use[_id] == 1:
            SeqIO.write(r, f_utr, 'fasta')
    for r in SeqIO.parse(open(cds2), 'fasta'):
        _id = r.id.split('|')[0]
        if to_use[_id] == 2:
            SeqIO.write(r, f_cds, 'fasta')
    for r in SeqIO.parse(open(pep2), 'fasta'):
        _id = r.id.split('|')[0]
        if to_use[_id] == 2:
            SeqIO.write(r, f_pep, 'fasta')
    for r in SeqIO.parse(open(utr2), 'fasta'):
        _id = r.id.split('|')[0]
        if to_use[_id] == 2:
            SeqIO.write(r, f_utr, 'fasta')

    f_cds.close()
    f_pep.close()
    f_utr.close()



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Program for picking best ORF from two ANGEL output files")
    parser.add_argument("input_prefix1", help="Input prefix for fileset 1 (ex: pabcio.ANGEL)")
    parser.add_argument("input_prefix2", help="Input prefix for fileset 2 (ex: genome_corrected.ANGEL")
    parser.add_argument("output_prefix", help="Output prefix")

    args = parser.parse_args()

    compare_files(args.input_prefix1, args.input_prefix2, args.output_prefix)