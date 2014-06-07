__author__ = 'lachesis'

import os, sys, subprocess

def sanity_check_cdhit():
    if os.system("cd-hit --version")!=256:
        raise Exception, "Cannot find cd-hit! Abort!"

def sanity_check_ATCG(func):
    def g(input_string, *arg):
        if type(input_string) is not str:
            raise TypeError, "Input string must be of string type!"
        for s in input_string:
            if s not in ('A','T','C','G'):
                raise TypeError, "Input string must consist of only ATCG!"
        return func(input_string, *arg)
    return g

@sanity_check_ATCG
def test(seq):
    print "just a test. input is", seq


def format_ORF_id(name, type, frame, start, end, strand):
    return "{name} type:{t} len:{l} strand:{strand} pos:{s}-{e}".format(\
        name=name, t=type, l=(end-start)/3, strand=strand, s=start+1, e=end)

def write_CDS_n_PEP(ORFs, output_prefix, min_utr_length=50, append_file=False, starting_index=1):
    """
    ORFs --- list of (Bio.SeqRecord, result, strand) where
             result is dict of frame --> list of (type, start, end)
    """
    index = starting_index
    f_cds = open(output_prefix + '.cds', 'w' if not append_file else 'a')
    f_pep = open(output_prefix + '.pep', 'w' if not append_file else 'a')
    f_utr = open(output_prefix + '.utr', 'w' if not append_file else 'a')
    for rec, result, strand in ORFs:
        seq_len = len(rec.seq)
        for frame, orfs in result.iteritems():
            for type, start, end in orfs:
                name = rec.id + '|m.' + str(index)
                index += 1
                if strand == '+':
                    orf_id = format_ORF_id(name, type, frame, start, end, strand)
                    f_cds.write(">{0}\n{1}\n".format(orf_id, rec.seq[start:end]))
                    f_pep.write(">{0}\n{1}\n".format(orf_id, rec.seq[start:end].translate()))
                    if start >= min_utr_length:
                        utr_id = format_ORF_id(name, '5UTR', 'NA', 0, start, strand)
                        f_utr.write(">{0}\n{1}\n".format(utr_id, rec.seq[:start]))
                    if seq_len - end + 1 >= min_utr_length:
                        utr_id = format_ORF_id(name, '3UTR', 'NA', end, seq_len, strand)
                        f_utr.write(">{0}\n{1}\n".format(utr_id, rec.seq[end:]))
                else: # strand == '-', need to adjust start, end, and seq
                    r_seq = rec.seq.reverse_complement()
                    orf_id = format_ORF_id(name, type, frame, seq_len-end, seq_len-start, strand)
                    f_cds.write(">{0}\n{1}\n".format(orf_id, r_seq[start:end]))
                    f_pep.write(">{0}\n{1}\n".format(orf_id, r_seq[start:end].translate()))
                    if start >= min_utr_length:  # 0:start
                        utr_id = format_ORF_id(name, '5UTR', 'NA', seq_len-start, seq_len-0, strand)
                        f_utr.write(">{0}\n{1}\n".format(utr_id, r_seq[:start]))
                    if seq_len - end + 1 >= min_utr_length:  #end:seq_len
                        utr_id = format_ORF_id(name, '3UTR', 'NA', 0, seq_len-end, strand)
                        f_utr.write(">{0}\n{1}\n".format(utr_id, r_seq[end:]))

    f_cds.close()
    f_pep.close()
    f_utr.close()



if __name__ == "__main__":
    test(sys.argv[1])