__author__ = 'etseng@pacificbiosciences.com'

import os, sys, subprocess
from Angel.ORFutils import sanity_check_ATCG, sanity_check_cdhit, write_CDS_n_PEP
import Angel.ORFscores as ORFscores
from Bio import SeqIO

@sanity_check_ATCG
def predict_longest_ORFs(seq, min_aa_length):
    """
    seq --- should be plain string in all upper case, A/T/C/G
    Return all longest ORFs that exceed <min_length>

    Returns: dict of <frame> --> list of (flag, <0-based start>, <1-based end>)
    NOTE that is the seq is reverse complemented, the handler function needs to rev the coords on its own
    """
    start_d, stop_d = ORFscores.find_start_stop_codons(seq)
    result = {0: [], 1: [], 2: []}

    n, m = len(seq)/3, len(seq)%3

    for frame in xrange(3):
        starts, stops = start_d[frame].keys(), stop_d[frame].keys()
        starts.sort()
        stops.sort()
        #print frame, starts, stops
        if len(stops) == 0: # no stop, so just output first (start, last)
            if len(starts) > 0 and n - starts[0] + 1 >= min_aa_length:
                result[frame].append(('dumb-3partial', starts[0]*3+frame, n*3+(frame if frame<=m else 0)))
        else: # has stop
            if len(starts) == 0: # 5' partial
                if  stops[0] + 1 >= min_aa_length:
                    result[frame].append(('dumb-5partial', frame, stops[0]*3+3+frame))
            else: # has at least one start and one stop
                i, j = 0, 0
                while j < len(stops):
                    if i == len(starts): break
                    if stops[j] - starts[i] + 1 >= min_aa_length:
                        result[frame].append(('dumb-complete', starts[i]*3+frame, stops[j]*3+3+frame))
                    j += 1 # move stop one step down
                    while i < len(starts) and starts[i] < stops[j-1]:
                        i += 1
                # check the very last possible ORF
                if i < len(starts) and (j == len(stops) or (j < len(stops) and starts[i] > stops[j])) and n - starts[i] + 1 >= min_aa_length:
                    result[frame].append(('dumb-3partial', starts[i]*3+frame, n*3+(frame if frame<=m else 0)))
    return result



def transdecoder_main(fasta_filename, output_prefix='dumb_orf', min_aa_length=100, use_rev_strand=False, use_top=500, cpus=8):
    sanity_check_cdhit()

    print >> sys.stderr, "predict longest ORFs...."
    # step 1. predict longest ORFs
    ORFs = [] # list of (sequence, result, strand)
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        seq = r.seq.tostring().upper()
        result = predict_longest_ORFs(seq, min_aa_length)
        ORFs.append((r, result, '+'))
        if use_rev_strand: # predict on - strand as well
            seq = r.seq.reverse_complement().tostring().upper()
            result = predict_longest_ORFs(seq, min_aa_length)
            ORFs.append((r, result, '-'))
    write_CDS_n_PEP(ORFs, output_prefix)

    print >> sys.stderr, "running CD-HIT to generate non-redundant set...."
    # step 2. use CD-hit to remove redundancy, then pick out top <use_top>
    cmd = "cd-hit -T {cpus} -i {o}.cds -o {o}.cds.nr90 -c 0.90 -n 5".format(o=output_prefix, cpus=cpus)
    subprocess.check_call(cmd, shell=True)

    lengths = [(len(r.seq), r) for r in SeqIO.parse(open(output_prefix+'.cds.nr90'), 'fasta')]
    lengths.sort(key=lambda x: x[0], reverse=True)
    lengths = lengths[:use_top]

    picked = []
    with open(output_prefix + '.training_' + str(use_top) + '.cds', 'w') as f:
        for _len, r in lengths:
            # ex: r.description >PB.1.1|chr1:26227060-26232896(-)|c242/f4p4/976|m.1 type:complete len:150 strand:+ pos:80-529
            f.write(">{0}\n{1}\n".format(r.description, r.seq))
            picked.append(r.id)

    with open(output_prefix + '.training_' + str(use_top) + '.utr', 'w') as f:
        for r in SeqIO.parse(open(output_prefix + '.utr'), 'fasta'):
            if r.id in picked:
                f.write(">{0}\n{1}\n".format(r.description, r.seq))






