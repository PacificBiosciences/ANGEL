__author__ = 'etseng@pacificbiosciences.com'

import os, sys, subprocess
import math, itertools
import random
from collections import defaultdict
from Angel.ORFutils import sanity_check_ATCG, sanity_check_cdhit, write_CDS_n_PEP, selective_write
import Angel.ORFscores as ORFscores
from Bio import SeqIO

@sanity_check_ATCG
def predict_longest_ORFs(seq, min_aa_length, use_firstORF=False):
    """
    seq --- should be plain string in all upper case, A/T/C/G
    Return the longest ORFs that exceed <min_length> (unless use_firstORF is True)

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
                result[frame].append(('dumb-3partial', starts[0]*3+frame, n*3+(frame if frame<=m else frame-3)))
        else: # has stop
            if len(starts) == 0: # 5' partial
                if  stops[0] + 1 >= min_aa_length:
                    result[frame].append(('dumb-5partial', frame, stops[0]*3+3+frame))
            else: # has at least one start and one stop
                i, j = 0, 0
                # if the first stop is smaller than i, find the first j s.t. stops[j-1] < start[0] < stops[j]
                if stops[0] < starts[0]:
                    while j < len(stops) and starts[0] < stops[j-1]:
                        j += 1
                # now: stops[j-1] < starts[0] < stops[j]
                while j < len(stops):
                    if i == len(starts): break
                    if stops[j] - starts[i] + 1 >= min_aa_length:
                        #rint frame, starts[i], stops[j]
                        result[frame].append(('dumb-complete', starts[i]*3+frame, stops[j]*3+3+frame))
                    j += 1 # move stop one step down
                    while i < len(starts) and starts[i] < stops[j-1]:
                        i += 1
                    # now starts[i] is between the last stop and this one
                # check the very last possible ORF
                if i < len(starts) and (j == len(stops) or (j < len(stops) and starts[i] > stops[j])) and n - starts[i] + 1 >= min_aa_length:
                    result[frame].append(('dumb-3partial', starts[i]*3+frame, n*3+(frame if frame<=m else frame-3)))

    # now pick the frame with the longest ORF!
    if all(len(v)==0 for v in result.itervalues()): # no ORF found
        return None


    best_frame, best_flag, best_s, best_e, best_len = None, None, None, None, 0
    if not use_firstORF: # find the longest ORF among all frames
        for _frame, v in result.iteritems():
            for (flag, s, e) in v:
                _len = e - s
                if _len > best_len:
                    best_frame, best_flag, best_s, best_e, best_len = \
                    _frame, flag, s, e, _len
    else: # use the first ORF among all frames
        for _frame, v in result.iteritems():
            for (flag, s, e) in v:
                _len = e - s
                if best_s is None or s < best_s or (s==best_s and _len>best_len):
                    best_frame, best_flag, best_s, best_e, best_len = \
                    _frame, flag, s, e, _len

    return {best_frame: [(best_flag, best_s, best_e)]}

def calculate_base_frequency(fasta_filename, output_filename, use_rev_strand=False):
    """
    Calculate A/T/G/C frequency based on <fasta_filename> and write to <output_filename> as:
    A <count> <freq> ...
    T
    C
    G
    Return: dict of nt --> frequency

    *NOTE: if rev strand if used, AT/CG will always be symmetrical
    """
    counts = {'A':0, 'T':0, 'C':0, 'G':0}
    freq = {}
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        seq = r.seq.tostring().upper()
        for s in seq: counts[s] += 1
        if use_rev_strand:
            seq = r.seq.reverse_complement().tostring().upper()
            for s in seq: counts[s] += 1
    _sum = sum(v for v in counts.itervalues())
    for k,v in counts.iteritems(): freq[k] = counts[k] * 1. / _sum
    with open(output_filename, 'w') as f:
        for k in ['A', 'T', 'C', 'G']:
            f.write("{0}\t{1}\t{2}\n".format(k, counts[k], freq[k]))
    return freq

def calculate_hexa_penta_score(cds_filename, base_freq, output_filename):
    """
    Return hexa, penta where each is
     frame --> hexamer or pentamer --> log score
    """
    log_scores = {} # key: (hexa-frame), value: log odds score
    hexamer = {0: defaultdict(lambda: 1), 1:defaultdict(lambda: 1), 2:defaultdict(lambda: 1)}
    pentamer = {0: defaultdict(lambda: 4), 1:defaultdict(lambda: 4), 2:defaultdict(lambda: 4)}
    for r in SeqIO.parse(open(cds_filename), 'fasta'):
        seq = r.seq.tostring().upper()
        seq_len = len(seq)
        #assert seq_len % 3 == 0
        for i in xrange(seq_len-5):
            frame = i % 3
            hexamer[frame][seq[i:i+6]] += 1
            pentamer[frame][seq[i:i+5]] += 1
        i = seq_len - 5
        pentamer[i%3][seq[i:i+5]] += 1

    f = open(output_filename, 'w')
    for hexa in itertools.product('ATCG', repeat=6):
        for frame in xrange(3):
            hexa = "".join(hexa)
            score = math.log(hexamer[frame][hexa]) - math.log(pentamer[frame][hexa[:5]]) - math.log(base_freq[hexa[-1]])
            f.write("{0}-{1}\t{2}\n".format(hexa, frame, score))
            log_scores[hexa+'-'+str(frame)] = score
    f.close()
    return log_scores

def score_cds_by_likelihood(cds_filename, log_scores):
    result = {} # seqid --> scores in all six frames
    f = open(cds_filename + '.scores', 'w')
    for r in SeqIO.parse(open(cds_filename), 'fasta'):
        scores = []
        seq = r.seq.tostring().upper()
        for frame in xrange(3):
            seq2 = seq[frame:]
            scores.append(sum(log_scores[seq2[i:i+6]+'-'+str(i%3)] for i in xrange(len(seq2)-5)))
        seq = r.seq.reverse_complement().tostring().upper()
        for frame in xrange(3):
            seq2 = seq[frame:]
            scores.append(sum(log_scores[seq2[i:i+6]+'-'+str(i%3)] for i in xrange(len(seq2)-5)))
        f.write("{0}\t{1}\n".format(r.id, "\t".join(map(str, scores))))
        result[r.id] = scores
    f.close()
    return result


def transdecoder_main(fasta_filename, output_prefix='dumb_orf', min_aa_length=100, use_rev_strand=False, use_firstORF=False, cpus=8):
    """
    1. Predict longest ORFs, write to <output_prefix>.cds|.utr|.pep
    2. Run CD-hit to get non-redundant set, then pick the top 500 for getting hexamer information, <output_prefix>.nr90.longest_500.cds
    3. Get base_freq out of <fasta_filename>, get hexamer scores out of (2)
    4. Score everything from (1) based on (3), write to <output_prefix>.cds.scores
    5. Output the FINAL <output_prefix>.final.cds|.utr|.pep based on the scores from (4)
    """
    sanity_check_cdhit()

    print >> sys.stderr, "predict longest ORFs...."
    # step 1. predict longest ORFs
    ORFs = [] # list of (sequence, result, strand)
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        seq = r.seq.tostring().upper()
        result = predict_longest_ORFs(seq, min_aa_length, use_firstORF) # result is {best_frame: [(best_flag, best_s, best_e)]}
        if result is not None:
            ORFs.append((r, result, '+'))
        if use_rev_strand: # predict on - strand as well
            seq = r.seq.reverse_complement().tostring().upper()
            result = predict_longest_ORFs(seq, min_aa_length, use_firstORF)
            if result is not None:
                ORFs.append((r, result, '-'))


    if use_firstORF: # no need to do scoring, just use firstORF
        # simply find the first ORF in ORFs
        write_CDS_n_PEP(ORFs, output_prefix + '.final')
        print >> sys.stderr, "Dumb ORF prediction done. Final output written to:", output_prefix + '.final.cds', \
            output_prefix + '.final.utr', output_prefix + '.final.pep'
        return # all done!
    else: # need to score, write this current one down first
        write_CDS_n_PEP(ORFs, output_prefix)

    print >> sys.stderr, "running CD-HIT to generate non-redundant set...."
    # step 2. use CD-hit to remove redundancy, then pick out top <use_top>
    cmd = "cd-hit -T {cpus} -M 0 -i {o}.cds -o {o}.nr90.cds -c 0.90 -n 5".format(o=output_prefix, cpus=cpus)
    subprocess.check_call(cmd, shell=True)

    lengths = [(len(r.seq), r) for r in SeqIO.parse(open(output_prefix+'.nr90.cds'), 'fasta')]
    lengths.sort(key=lambda x: x[0], reverse=True)
    lengths = lengths[:500]

    cds_nr_selected_filename = output_prefix + '.nr90.longest_500.cds'
    with open(cds_nr_selected_filename, 'w') as f:
        for _len, r in lengths:
            # ex: r.description >PB.1.1|chr1:26227060-26232896(-)|c242/f4p4/976|m.1 type:complete len:150 strand:+ pos:80-529
            f.write(">{0}\n{1}\n".format(r.description, r.seq))
    print >> sys.stderr, "Longest 500 non-redundant predicted ORFs written to:", cds_nr_selected_filename


    # step 3. get base_freq & hexamer scores
    print >> sys.stderr, "Calculating base frequency from", fasta_filename
    base_freq = calculate_base_frequency(fasta_filename, fasta_filename+'.base_freq', use_rev_strand)
    print >> sys.stderr, "Calculating hexamer scores from", cds_nr_selected_filename
    log_scores = calculate_hexa_penta_score(cds_nr_selected_filename, base_freq, cds_nr_selected_filename+'.hexamer.scores')

    # step 4. score all predicted longest ORFs using log score
    print >> sys.stderr, "Scoring predicted ORFs...."
    scored_result = score_cds_by_likelihood(output_prefix + '.cds', log_scores)

    # step 5. output FINAL, where longest ORFs are output ONLY if its score in frame0 is higher than all other 5
    picked_ids = []
    for rec_seq_id, scores in scored_result.iteritems():
        if scores[0] > 0 and scores[0] == max(scores):
            picked_ids.append(rec_seq_id)

    selective_write(output_prefix + '.cds', output_prefix + '.final.cds', picked_ids)
    selective_write(output_prefix + '.utr', output_prefix + '.final.utr', picked_ids)
    selective_write(output_prefix + '.pep', output_prefix + '.final.pep', picked_ids)

    print >> sys.stderr, "Dumb ORF prediction done. Final output written to:", output_prefix + '.final.cds', \
        output_prefix + '.final.utr', output_prefix + '.final.pep'


def select_for_training(input_prefix, output_prefix, use_top=500, choose_random=False, cpus=8):
    """
    <input_prefix>.cds|.pep|.utr, must exist! Probably the output from transdecoder_main()

    use_top --- number of top records to use
    random --- if True, choose randomly instead of the longest top ones
    """
    print >> sys.stderr, "running CD-HIT to generate non-redundant set...."
    cmd = "cd-hit -T {cpus} -M 0 -i {o}.cds -o {o}.nr90.cds -c 0.90 -n 5".format(o=input_prefix, cpus=cpus)
    subprocess.check_call(cmd, shell=True)

    lengths = [(len(r.seq), r.id) for r in SeqIO.parse(open(input_prefix+'.nr90.cds'), 'fasta')]
    if not choose_random:
        print >> sys.stderr, "Selecting longest {0} entries from non-redundant set....".format(use_top)
        lengths.sort(key=lambda x: x[0], reverse=True)
        lengths = lengths[:use_top]
    else:
        print >> sys.stderr, "Selecting random {0} entries from non-redundant set....".format(use_top)
        lengths = random.sample(lengths, min(len(lengths), use_top))

    picked_ids = [ rec_seq_id for (seq_len, rec_seq_id) in lengths ]
    selective_write(input_prefix + '.cds', output_prefix + '.cds', picked_ids)
    selective_write(input_prefix + '.utr', output_prefix + '.utr', picked_ids)
    selective_write(input_prefix + '.pep', output_prefix + '.pep', picked_ids)