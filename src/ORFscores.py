__author__ = 'etseng@pacificbiosciences.com'

import os, sys
import itertools
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import Bio.Data.CodonTable as CodonTable
import Angel.c_ORFscores as sp2
from multiprocessing import Process
import Angel.findPath as findPath

AMINO_LETTERS = 'ACDEFGHIKLMNPQRSTVWY*'
NT_LETTERS = 'GATC'
CODON_LETTERS = ["GGG","GGA","GGT","GGC","GAG","GAA","GAT","GAC","GTG","GTA","GTT","GTC","GCG","GCA","GCT","GCC","AGG","AGA","AGT","AGC","AAG","AAA","AAT","AAC","ATG","ATA","ATT","ATC","ACG","ACA","ACT","ACC","TGG","TGA","TGT","TGC","TAG","TAA","TAT","TAC","TTG","TTA","TTT","TTC","TCG","TCA","TCT","TCC","CGG","CGA","CGT","CGC","CAG","CAA","CAT","CAC","CTG","CTA","CTT","CTC","CCG","CCA","CCT","CCC"]
DIAMINO_LETTERS = [\
    ('A', 'A'),('A', 'C'),('A', 'D'),('A', 'E'),('A', 'F'),('A', 'G'),('A', 'H'),('A', 'I'),('A', 'K'),('A', 'L'),('A', 'M'),('A', 'N'),('A', 'P'),('A', 'Q'),('A', 'R'),('A', 'S'),('A', 'T'),('A', 'V'),('A', 'W'),('A', 'Y'),('A', '*'),('C', 'A'),('C', 'C'),('C', 'D'),('C', 'E'),('C', 'F'),('C', 'G'),('C', 'H'),('C', 'I'),('C', 'K'),('C', 'L'),('C', 'M'),('C', 'N'),('C', 'P'),('C', 'Q'),('C', 'R'),('C', 'S'),('C', 'T'),('C', 'V'),('C', 'W'),('C', 'Y'),('C', '*'),('D', 'A'),('D', 'C'),('D', 'D'),('D', 'E'),('D', 'F'),('D', 'G'),('D', 'H'),('D', 'I'),('D', 'K'),('D', 'L'),('D', 'M'),('D', 'N'),('D', 'P'),('D', 'Q'),('D', 'R'),('D', 'S'),('D', 'T'),('D', 'V'),('D', 'W'),('D', 'Y'),('D', '*'),('E', 'A'),('E', 'C'),('E', 'D'),('E', 'E'),('E', 'F'),('E', 'G'),('E', 'H'),('E', 'I'),('E', 'K'),('E', 'L'),('E', 'M'),('E', 'N'),('E', 'P'),('E', 'Q'),('E', 'R'),('E', 'S'),('E', 'T'),('E', 'V'),('E', 'W'),('E', 'Y'),('E', '*'),('F', 'A'),('F', 'C'),('F', 'D'),('F', 'E'),('F', 'F'),('F', 'G'),('F', 'H'),\
    ('F', 'I'),('F', 'K'),('F', 'L'),('F', 'M'),('F', 'N'),('F', 'P'),('F', 'Q'),('F', 'R'),('F', 'S'),('F', 'T'),('F', 'V'),('F', 'W'),('F', 'Y'),('F', '*'),('G', 'A'),('G', 'C'),('G', 'D'),('G', 'E'),('G', 'F'),('G', 'G'),('G', 'H'),('G', 'I'),('G', 'K'),('G', 'L'),('G', 'M'),('G', 'N'),('G', 'P'),('G', 'Q'),('G', 'R'),('G', 'S'),('G', 'T'),('G', 'V'),('G', 'W'),('G', 'Y'),('G', '*'),('H', 'A'),('H', 'C'),('H', 'D'),('H', 'E'),('H', 'F'),('H', 'G'),('H', 'H'),('H', 'I'),('H', 'K'),('H', 'L'),('H', 'M'),('H', 'N'),('H', 'P'),('H', 'Q'),('H', 'R'),('H', 'S'),('H', 'T'),('H', 'V'),('H', 'W'),('H', 'Y'),('H', '*'),('I', 'A'),('I', 'C'),('I', 'D'),('I', 'E'),('I', 'F'),('I', 'G'),('I', 'H'),('I', 'I'),('I', 'K'),('I', 'L'),('I', 'M'),('I', 'N'),('I', 'P'),('I', 'Q'),('I', 'R'),('I', 'S'),('I', 'T'),('I', 'V'),('I', 'W'),('I', 'Y'),('I', '*'),('K', 'A'),('K', 'C'),('K', 'D'),('K', 'E'),('K', 'F'),('K', 'G'),('K', 'H'),('K', 'I'),('K', 'K'),('K', 'L'),('K', 'M'),('K', 'N'),('K', 'P'),('K', 'Q'),\
    ('K', 'R'),('K', 'S'),('K', 'T'),('K', 'V'),('K', 'W'),('K', 'Y'),('K', '*'),('L', 'A'),('L', 'C'),('L', 'D'),('L', 'E'),('L', 'F'),('L', 'G'),('L', 'H'),('L', 'I'),('L', 'K'),('L', 'L'),('L', 'M'),('L', 'N'),('L', 'P'),('L', 'Q'),('L', 'R'),('L', 'S'),('L', 'T'),('L', 'V'),('L', 'W'),('L', 'Y'),('L', '*'),('M', 'A'),('M', 'C'),('M', 'D'),('M', 'E'),('M', 'F'),('M', 'G'),('M', 'H'),('M', 'I'),('M', 'K'),('M', 'L'),('M', 'M'),('M', 'N'),('M', 'P'),('M', 'Q'),('M', 'R'),('M', 'S'),('M', 'T'),('M', 'V'),('M', 'W'),('M', 'Y'),('M', '*'),('N', 'A'),('N', 'C'),('N', 'D'),('N', 'E'),('N', 'F'),('N', 'G'),('N', 'H'),('N', 'I'),('N', 'K'),('N', 'L'),('N', 'M'),('N', 'N'),('N', 'P'),('N', 'Q'),('N', 'R'),('N', 'S'),('N', 'T'),('N', 'V'),('N', 'W'),('N', 'Y'),('N', '*'),('P', 'A'),('P', 'C'),('P', 'D'),('P', 'E'),('P', 'F'),('P', 'G'),('P', 'H'),('P', 'I'),('P', 'K'),('P', 'L'),('P', 'M'),('P', 'N'),('P', 'P'),('P', 'Q'),('P', 'R'),('P', 'S'),('P', 'T'),('P', 'V'),('P', 'W'),('P', 'Y'),('P', '*'),('Q', 'A'),('Q', 'C'),('Q', 'D'),('Q', 'E'),('Q', 'F'),('Q', 'G'),('Q', 'H'),('Q', 'I'),('Q', 'K'),('Q', 'L'),('Q', 'M'),('Q', 'N'),('Q', 'P'),('Q', 'Q'),('Q', 'R'),('Q', 'S'),('Q', 'T'),('Q', 'V'),('Q', 'W'),('Q', 'Y'),('Q', '*'),('R', 'A'),('R', 'C'),('R', 'D'),('R', 'E'),('R', 'F'),('R', 'G'),('R', 'H'),('R', 'I'),('R', 'K'),('R', 'L'),('R', 'M'),('R', 'N'),('R', 'P'),('R', 'Q'),('R', 'R'),('R', 'S'),('R', 'T'),('R', 'V'),('R', 'W'),('R', 'Y'),('R', '*'),('S', 'A'),('S', 'C'),('S', 'D'),('S', 'E'),('S', 'F'),('S', 'G'),('S', 'H'),('S', 'I'),('S', 'K'),('S', 'L'),('S', 'M'),('S', 'N'),('S', 'P'),('S', 'Q'),('S', 'R'),('S', 'S'),('S', 'T'),('S', 'V'),('S', 'W'),('S', 'Y'),('S', '*'),('T', 'A'),('T', 'C'),('T', 'D'),('T', 'E'),('T', 'F'),('T', 'G'),('T', 'H'),('T', 'I'),('T', 'K'),('T', 'L'),('T', 'M'),('T', 'N'),('T', 'P'),('T', 'Q'),('T', 'R'),('T', 'S'),('T', 'T'),('T', 'V'),('T', 'W'),('T', 'Y'),('T', '*'),('V', 'A'),('V', 'C'),('V', 'D'),('V', 'E'),('V', 'F'),('V', 'G'),('V', 'H'),('V', 'I'),('V', 'K'),('V', 'L'),('V', 'M'),('V', 'N'),('V', 'P'),('V', 'Q'),('V', 'R'),('V', 'S'),('V', 'T'),('V', 'V'),('V', 'W'),('V', 'Y'),('V', '*'),('W', 'A'),('W', 'C'),('W', 'D'),('W', 'E'),('W', 'F'),('W', 'G'),('W', 'H'),('W', 'I'),('W', 'K'),('W', 'L'),('W', 'M'),('W', 'N'),('W', 'P'),('W', 'Q'),('W', 'R'),('W', 'S'),('W', 'T'),('W', 'V'),('W', 'W'),('W', 'Y'),('W', '*'),('Y', 'A'),('Y', 'C'),('Y', 'D'),('Y', 'E'),('Y', 'F'),('Y', 'G'),('Y', 'H'),('Y', 'I'),('Y', 'K'),('Y', 'L'),('Y', 'M'),('Y', 'N'),('Y', 'P'),('Y', 'Q'),('Y', 'R'),('Y', 'S'),('Y', 'T'),('Y', 'V'),('Y', 'W'),('Y', 'Y'),('Y', '*'),('*', 'A'),('*', 'C'),('*', 'D'),('*', 'E'),('*', 'F'),('*', 'G'),('*', 'H'),('*', 'I'),('*', 'K'),('*', 'L'),('*', 'M'),('*', 'N'),('*', 'P'),('*', 'Q'),('*', 'R'),('*', 'S'),('*', 'T'),('*', 'V'),('*', 'W'),('*', 'Y'),('*', '*')]


START_CODONS = ['ATG']  # for now, stick with most canonical start codon
STOP_CODONS = ['TAA', 'TAG', 'TGA']


def make_amino_scores(freq):
    arr = []
    # make sure it's in the right *order*
    for x in AMINO_LETTERS:
        arr.append(freq[x])
    return arr

def make_diamino_scores(freq, diamino_range):
    """
    In the order of (k, A_i, B_j)
    """
    arr = [None] * (len(diamino_range) * len(DIAMINO_LETTERS))
    i = 0
    for k in diamino_range:
        for x,y in DIAMINO_LETTERS:            
            arr[i] = freq[k][(x, y)]
            i += 1
    return arr

def make_codon_scores(aa_freq, codon_freq):
    arr = []
    letters = NT_LETTERS
    for codon, aa in CodonTable.standard_dna_table.forward_table.iteritems():
        arr.append(aa_freq[aa]*1./codon_freq[codon])
    for codon in CodonTable.standard_dna_table.stop_codons:
        arr.append(aa_freq['*']*1./codon_freq[codon])
    return arr


def make_data_smart(seq, pseudo, window_size=96, step_size=3, frame_shift=0):
    seq = seq.upper()
    ss = str(seq[frame_shift:])
    a, b = len(ss)/3, len(ss)%3
    aa_seq_end = a * 3 + frame_shift if frame_shift <= b else frame_shift-3
    aa = str(seq[frame_shift:aa_seq_end].translate())
    aa_window_size = window_size / 3 
    n = len(ss)
    
    o = sp2.CDSWindowFeat()
    cur_a = aa[:aa_window_size]
    cur_s = ss[:window_size]
    o.calc_amino_count(cur_a)
    o.calc_diamino_count(cur_a, len(cur_a))
    o.calc_codon_count(cur_s, len(cur_s))
    aa_freq = o.get_amino_freq(pseudo, .0001)
    di_freq = o.get_diamino_freq(pseudo, .0001)
    codon_freq = o.get_codon_freq(pseudo, .0001)
    arr = make_amino_scores(aa_freq)+make_diamino_scores(di_freq,o.diamino_range)+make_codon_scores(aa_freq,codon_freq)
    data = [arr]    
    for i in xrange(1, len(aa) - aa_window_size):
        # s advances by 3, now at ss[i*3:i*3+window_size]
        # a advances by 1, now at aa[i:i+aa_window_size]
        o.calc_amino_count(cur_a[0], -1)
        o.calc_amino_count(aa[i+aa_window_size-1], 1)
        o.calc_codon_count(cur_s[:3], 3, -1)
        o.calc_codon_count(ss[(i-1)*3+window_size:(i)*3+window_size], 3, 1)
        o.diamino_changed = {} # must clear out THIS! 
        o.deduct_diamino_count(cur_a, len(cur_a), 1)        
        cur_a = aa[i:i+aa_window_size]        
        o.add_diamino_count(cur_a, len(cur_a), 1)
        cur_s = ss[i*3:i*3+window_size]
        aa_freq = o.get_amino_freq(pseudo, .0001)
        di_freq = o.get_diamino_freq(pseudo, .0001, False)
        codon_freq = o.get_codon_freq(pseudo, .0001)
        arr = make_amino_scores(aa_freq)+make_diamino_scores(di_freq,o.diamino_range)+make_codon_scores(aa_freq,codon_freq)
        data.append(arr)
        #if i > 5: break
    return data

def find_start_stop_codons(seq):
    """
    Given a sequence, return a dict of:
      frame (0, 1, 2) --> n-th codon --> "TAA/TAG/TGA"
      frame (0, 1, 2) --> n-th codon --> "TAG"

    Returns: start_dict, stop_dict
    """
    stop_d = {0: {}, 1: {}, 2: {}}
    start_d = {0: {}, 1: {}, 2: {}}
    for i in xrange(len(seq)-2):
        if seq[i:i+3] in STOP_CODONS:
            frame = i % 3
            nth = (i-frame) / 3
            stop_d[i%3][nth] = seq[i:i+3]
        elif seq[i:i+3] in START_CODONS:
            frame = i % 3
            nth = (i-frame) / 3
            start_d[i%3][nth] = seq[i:i+3]
    return start_d, stop_d


# def find_indel(seq, min_i, max_i, bdt, background, target_frame, old_score, stop_when_found=True):
#     best_score, best_moves = old_score, []
#     # try deleting one base first
#     for i in xrange(min_i*3, max_i*3):
#         new_seq = seq[:i] + seq[(i+1):]
#         start_dict, stop_dict = find_start_stop_codons(new_seq.tostring())
#         if len(stop_dict[target_frame]) > 0: continue # not feasible if introduces a stop codon
#         stuff = make_data(new_seq, background, frame_shift=target_frame)
#         score = sum(bdt.predict(stuff))
#         if score > best_score:
#             best_score, best_moves = score, [(i, new_seq)]
#             if stop_when_found: return best_score, best_moves
#         elif score == best_score:
#             best_moves.append((i, new_seq))
#
#     for i in xrange(min_i*3, max_i*3):
#         new_seq = seq[:i] + seq[(i+2):]
#         start_dict, stop_dict = find_start_stop_codons(new_seq.tostring())
#         if len(stop_dict[target_frame]) > 0: continue # not feasible if introduces a stop codon
#         stuff = make_data(new_seq, background, frame_shift=target_frame)
#         score = sum(bdt.predict(stuff))
#         if score > best_score:
#             best_score, best_moves = score, [((i,i+1), new_seq)]
#             if stop_when_found: return best_score, best_moves
#         elif score == best_score:
#             best_moves.append(((i,i+1), new_seq))
#
#     return best_score, best_moves

def find_chunks(rec, bdt, o_all):
    """
    rec --- Bio.SeqRecord
    bdt --- predictor object
    o_all --- background freq
    """
    rec.seq = rec.seq.upper() # need to make sure all capitalized
    stuff0 = make_data_smart(rec.seq, o_all, frame_shift=0)
    stuff1 = make_data_smart(rec.seq, o_all, frame_shift=1)
    stuff2 = make_data_smart(rec.seq, o_all, frame_shift=2)
    ans0 = bdt.predict(stuff0)
    ans1 = bdt.predict(stuff1)
    ans2 = bdt.predict(stuff2)
    start_dict,stop_dict = find_start_stop_codons(str(rec.seq))
    E,O = findPath.make_DP_matrix(ans0, ans1, ans2, start_dict, stop_dict)

    i, j = E.argmax()/3, E.argmax()%3
    last = i, j
    chunks = [] # (frame, start, stop)
    while True:
        try:
            i, j = O[(i, j)]
        except:
            chunks.append((last[1], last[0], i))
            break
        if j!=last[1]:
            #print "{0}-{1}, frame{2}".format(last[0],i, last[1])
            chunks.append((last[1], last[0], i))
            last = i, j
    return chunks, E, O, ans0, ans1, ans2, start_dict, stop_dict




def extend_in_frame(cur_begin, cur_end, starts, stops, size):
    """
    starts --- should be start_dict[frame], which is position --> "ATG" etc
    stops --- should be stop_dict[frame], which is position --> "TAG" etc

    Find the furthest valid begin, which is either the beginning of the seq (incomplete 5')
    or the furthest start codon, with no stop codon in between
    Similarity, find closest valid end, which is either end of the seq (incomplete 3', unlikely)
    or the closest stop codon
    """
    i = cur_begin + 1
    furthest_start = None
    j = cur_end - 1
    closest_stop = None
    while i >= 0:
        if i in starts:
            furthest_start = i
        if i in stops:
            furthest_start = i + 1
            # find instead the furthest possible stop, going in increment
            for k in xrange(i+1, cur_end):
                if k in starts:
                    furthest_start = k
                    break
            break
        i -= 1
    while j <= size:
        if j in stops:
            closest_stop = j
            break
        j += 1

    return furthest_start, closest_stop

def predict_ORF(rec, bdt, o_all, min_aa_len=200):
    """
    (1) if only a single chunk predicted --- confident ORF
    (2) if more than one chunk but only one with > min size, and proper stop --- semi-confident
    (3) if more than one chunk but with multiple > min size --- suspicious
    """
    chunks, E, O, ans0, ans1, ans2, start_dict, stop_dict = find_chunks(rec, bdt, o_all)
    ans = {0: ans0, 1: ans1, 2: ans2}

    n, m = len(rec.seq)/3, len(rec.seq)%3
    if len(chunks) == 1:
        _frame, _end, _begin = chunks[0]
        e_start, e_stop = extend_in_frame(_begin, _end, start_dict[_frame], stop_dict[_frame], (len(rec.seq)-_frame)/3)
        #print "single chunk for", rec.id
        #print chunks[0], "extended to", e_start, e_stop

        if e_stop is None:
            if e_start is None:
                flag = "confident-internal"
                aa_len = n + 1
            else:
                flag = "confident-3partial"
                aa_len = n - e_start + 1
            #return "confident-3partial", rec.id, [(_frame, e_stop, e_start)]
        elif e_start is None or e_start not in start_dict[_frame]:
            flag = "confident-5partial"
            aa_len = e_stop + 1 if e_start is None else (e_stop-e_start+1)
            #return 'confident-5partial', rec.id, [(_frame, e_stop, e_start)]
        else:
            flag = "confident-complete"
            aa_len = e_stop - e_start + 1
            #return "confident-complete", rec.id, [(_frame, e_stop, e_start)]
        if aa_len >= min_aa_len:
            return flag, rec.id, [(_frame, e_stop, e_start)]
    else:
        good = []
        for _frame,_end,_begin in chunks:
            e_start, e_stop = extend_in_frame(_begin, _end, start_dict[_frame], stop_dict[_frame], (len(rec.seq)-_frame)/3)
            e_start2 = 0 if e_start is None else e_start
            e_stop2 = (len(rec.seq)-_frame) / 3 if e_stop is None else e_stop
            e_len = e_stop2 - e_start2 + 1
            if e_len >= min_aa_len and sum(ans[_frame][e_start2:e_stop2]) >= .5 * e_len:
                good.append((_frame, e_stop, e_start))
        if len(good) == 1:
            #print "single good for", rec.id, good[0]
            return "likely-NA", rec.id, good
        else:
            #print "multi chunk"
            return "suspicious-NA", rec.id, good

    return None, None, []
    #return chunks, E, O, ans0, ans1, ans2, start_dict, stop_dict


def main_for_Gloria(input_fasta, output_prefix, bdt, o_all):
    f_csv = open(output_prefix+'.csv', 'w')
    f_csv.write("id\tstatus\tcompleteness\tframe\trange\n")
    f_fa = open(output_prefix+'.fa', 'w')
    f_pep = open(output_prefix+'.pep', 'w')
    for rec in SeqIO.parse(open(input_fasta), 'fasta'):
        print >> sys.stderr, "predicting for", rec.id
        flag, name, good = predict_ORF(rec, bdt, o_all, min_chunk_size=50)
        #print >> sys.stderr, flag, name, good
        status, completeness = flag.split('-')
        for _frame, _stop, _start in good:
            s = _start * 3 + _frame if _start is not None else _frame
            e = _stop * 3 + _frame + 3 if _stop is not None else (len(rec.seq)/3)*3 + (_frame if len(rec.seq)%3 >= _frame else 0)
            f_csv.write("{id}\t{stat}\t{com}\t{fr}\t{s}-{e}\n".format(\
                id=rec.id, stat=status, com=completeness, fr=_frame, s=s+1, e=e))
            f_fa.write(">{id};{stat};{com};{fr};{s}-{e}\n".format(\
                id=rec.id, stat=status, com=completeness, fr=_frame, s=s+1, e=e))
            f_fa.write("{seq}\n".format(seq=rec.seq[s:e]))
            f_pep.write(">{id};{stat};{com};{fr};{s}-{e}\n".format(\
                id=rec.id, stat=status, com=completeness, fr=_frame, s=s+1, e=e))
            f_pep.write("{seq}\n".format(seq=rec.seq[s:e].translate()))
    f_csv.close()
    f_fa.close()
    f_pep.close()

def distribute_forGloria(list_of_fasta, bdt, o_all):
    workers = []
    for input_fasta in list_of_fasta:
        print >> sys.stderr, "Pool worker for", input_fasta
        p = Process(target=main_for_Gloria, args=(input_fasta, input_fasta+'.ORF', bdt, o_all,))
        p.start()
        workers.append(p)

    for p in workers:
        print >> sys.stderr, "waiting for worker", p.name
        p.join()
