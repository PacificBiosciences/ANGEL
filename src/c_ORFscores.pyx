__author__ = 'etseng@pacificbiosciences.com'

import itertools
from Bio.Alphabet import IUPAC
import Bio.Data.CodonTable as CodonTable

from libcpp.utility cimport pair
from libcpp.map cimport map 
from libcpp.string cimport string
from libcpp cimport bool

AMINO_LETTERS = 'ACDEFGHIKLMNPQRSTVWY*'
NT_LETTERS = 'GATC'
CODON_LETTERS = ["GGG","GGA","GGT","GGC","GAG","GAA","GAT","GAC","GTG","GTA","GTT","GTC","GCG","GCA","GCT","GCC","AGG","AGA","AGT","AGC","AAG","AAA","AAT","AAC","ATG","ATA","ATT","ATC","ACG","ACA","ACT","ACC","TGG","TGA","TGT","TGC","TAG","TAA","TAT","TAC","TTG","TTA","TTT","TTC","TCG","TCA","TCT","TCC","CGG","CGA","CGT","CGC","CAG","CAA","CAT","CAC","CTG","CTA","CTT","CTC","CCG","CCA","CCT","CCC"]
DIAMINO_LETTERS = [\
    ('A', 'A'),('A', 'C'),('A', 'D'),('A', 'E'),('A', 'F'),('A', 'G'),('A', 'H'),('A', 'I'),('A', 'K'),('A', 'L'),('A', 'M'),('A', 'N'),('A', 'P'),('A', 'Q'),('A', 'R'),('A', 'S'),('A', 'T'),('A', 'V'),('A', 'W'),('A', 'Y'),('A', '*'),('C', 'A'),('C', 'C'),('C', 'D'),('C', 'E'),('C', 'F'),('C', 'G'),('C', 'H'),('C', 'I'),('C', 'K'),('C', 'L'),('C', 'M'),('C', 'N'),('C', 'P'),('C', 'Q'),('C', 'R'),('C', 'S'),('C', 'T'),('C', 'V'),('C', 'W'),('C', 'Y'),('C', '*'),('D', 'A'),('D', 'C'),('D', 'D'),('D', 'E'),('D', 'F'),('D', 'G'),('D', 'H'),('D', 'I'),('D', 'K'),('D', 'L'),('D', 'M'),('D', 'N'),('D', 'P'),('D', 'Q'),('D', 'R'),('D', 'S'),('D', 'T'),('D', 'V'),('D', 'W'),('D', 'Y'),('D', '*'),('E', 'A'),('E', 'C'),('E', 'D'),('E', 'E'),('E', 'F'),('E', 'G'),('E', 'H'),('E', 'I'),('E', 'K'),('E', 'L'),('E', 'M'),('E', 'N'),('E', 'P'),('E', 'Q'),('E', 'R'),('E', 'S'),('E', 'T'),('E', 'V'),('E', 'W'),('E', 'Y'),('E', '*'),('F', 'A'),('F', 'C'),('F', 'D'),('F', 'E'),('F', 'F'),('F', 'G'),('F', 'H'),\
    ('F', 'I'),('F', 'K'),('F', 'L'),('F', 'M'),('F', 'N'),('F', 'P'),('F', 'Q'),('F', 'R'),('F', 'S'),('F', 'T'),('F', 'V'),('F', 'W'),('F', 'Y'),('F', '*'),('G', 'A'),('G', 'C'),('G', 'D'),('G', 'E'),('G', 'F'),('G', 'G'),('G', 'H'),('G', 'I'),('G', 'K'),('G', 'L'),('G', 'M'),('G', 'N'),('G', 'P'),('G', 'Q'),('G', 'R'),('G', 'S'),('G', 'T'),('G', 'V'),('G', 'W'),('G', 'Y'),('G', '*'),('H', 'A'),('H', 'C'),('H', 'D'),('H', 'E'),('H', 'F'),('H', 'G'),('H', 'H'),('H', 'I'),('H', 'K'),('H', 'L'),('H', 'M'),('H', 'N'),('H', 'P'),('H', 'Q'),('H', 'R'),('H', 'S'),('H', 'T'),('H', 'V'),('H', 'W'),('H', 'Y'),('H', '*'),('I', 'A'),('I', 'C'),('I', 'D'),('I', 'E'),('I', 'F'),('I', 'G'),('I', 'H'),('I', 'I'),('I', 'K'),('I', 'L'),('I', 'M'),('I', 'N'),('I', 'P'),('I', 'Q'),('I', 'R'),('I', 'S'),('I', 'T'),('I', 'V'),('I', 'W'),('I', 'Y'),('I', '*'),('K', 'A'),('K', 'C'),('K', 'D'),('K', 'E'),('K', 'F'),('K', 'G'),('K', 'H'),('K', 'I'),('K', 'K'),('K', 'L'),('K', 'M'),('K', 'N'),('K', 'P'),('K', 'Q'),\
    ('K', 'R'),('K', 'S'),('K', 'T'),('K', 'V'),('K', 'W'),('K', 'Y'),('K', '*'),('L', 'A'),('L', 'C'),('L', 'D'),('L', 'E'),('L', 'F'),('L', 'G'),('L', 'H'),('L', 'I'),('L', 'K'),('L', 'L'),('L', 'M'),('L', 'N'),('L', 'P'),('L', 'Q'),('L', 'R'),('L', 'S'),('L', 'T'),('L', 'V'),('L', 'W'),('L', 'Y'),('L', '*'),('M', 'A'),('M', 'C'),('M', 'D'),('M', 'E'),('M', 'F'),('M', 'G'),('M', 'H'),('M', 'I'),('M', 'K'),('M', 'L'),('M', 'M'),('M', 'N'),('M', 'P'),('M', 'Q'),('M', 'R'),('M', 'S'),('M', 'T'),('M', 'V'),('M', 'W'),('M', 'Y'),('M', '*'),('N', 'A'),('N', 'C'),('N', 'D'),('N', 'E'),('N', 'F'),('N', 'G'),('N', 'H'),('N', 'I'),('N', 'K'),('N', 'L'),('N', 'M'),('N', 'N'),('N', 'P'),('N', 'Q'),('N', 'R'),('N', 'S'),('N', 'T'),('N', 'V'),('N', 'W'),('N', 'Y'),('N', '*'),('P', 'A'),('P', 'C'),('P', 'D'),('P', 'E'),('P', 'F'),('P', 'G'),('P', 'H'),('P', 'I'),('P', 'K'),('P', 'L'),('P', 'M'),('P', 'N'),('P', 'P'),('P', 'Q'),('P', 'R'),('P', 'S'),('P', 'T'),('P', 'V'),('P', 'W'),('P', 'Y'),('P', '*'),('Q', 'A'),('Q', 'C'),('Q', 'D'),('Q', 'E'),('Q', 'F'),('Q', 'G'),('Q', 'H'),('Q', 'I'),('Q', 'K'),('Q', 'L'),('Q', 'M'),('Q', 'N'),('Q', 'P'),('Q', 'Q'),('Q', 'R'),('Q', 'S'),('Q', 'T'),('Q', 'V'),('Q', 'W'),('Q', 'Y'),('Q', '*'),('R', 'A'),('R', 'C'),('R', 'D'),('R', 'E'),('R', 'F'),('R', 'G'),('R', 'H'),('R', 'I'),('R', 'K'),('R', 'L'),('R', 'M'),('R', 'N'),('R', 'P'),('R', 'Q'),('R', 'R'),('R', 'S'),('R', 'T'),('R', 'V'),('R', 'W'),('R', 'Y'),('R', '*'),('S', 'A'),('S', 'C'),('S', 'D'),('S', 'E'),('S', 'F'),('S', 'G'),('S', 'H'),('S', 'I'),('S', 'K'),('S', 'L'),('S', 'M'),('S', 'N'),('S', 'P'),('S', 'Q'),('S', 'R'),('S', 'S'),('S', 'T'),('S', 'V'),('S', 'W'),('S', 'Y'),('S', '*'),('T', 'A'),('T', 'C'),('T', 'D'),('T', 'E'),('T', 'F'),('T', 'G'),('T', 'H'),('T', 'I'),('T', 'K'),('T', 'L'),('T', 'M'),('T', 'N'),('T', 'P'),('T', 'Q'),('T', 'R'),('T', 'S'),('T', 'T'),('T', 'V'),('T', 'W'),('T', 'Y'),('T', '*'),('V', 'A'),('V', 'C'),('V', 'D'),('V', 'E'),('V', 'F'),('V', 'G'),('V', 'H'),('V', 'I'),('V', 'K'),('V', 'L'),('V', 'M'),('V', 'N'),('V', 'P'),('V', 'Q'),('V', 'R'),('V', 'S'),('V', 'T'),('V', 'V'),('V', 'W'),('V', 'Y'),('V', '*'),('W', 'A'),('W', 'C'),('W', 'D'),('W', 'E'),('W', 'F'),('W', 'G'),('W', 'H'),('W', 'I'),('W', 'K'),('W', 'L'),('W', 'M'),('W', 'N'),('W', 'P'),('W', 'Q'),('W', 'R'),('W', 'S'),('W', 'T'),('W', 'V'),('W', 'W'),('W', 'Y'),('W', '*'),('Y', 'A'),('Y', 'C'),('Y', 'D'),('Y', 'E'),('Y', 'F'),('Y', 'G'),('Y', 'H'),('Y', 'I'),('Y', 'K'),('Y', 'L'),('Y', 'M'),('Y', 'N'),('Y', 'P'),('Y', 'Q'),('Y', 'R'),('Y', 'S'),('Y', 'T'),('Y', 'V'),('Y', 'W'),('Y', 'Y'),('Y', '*'),('*', 'A'),('*', 'C'),('*', 'D'),('*', 'E'),('*', 'F'),('*', 'G'),('*', 'H'),('*', 'I'),('*', 'K'),('*', 'L'),('*', 'M'),('*', 'N'),('*', 'P'),('*', 'Q'),('*', 'R'),('*', 'S'),('*', 'T'),('*', 'V'),('*', 'W'),('*', 'Y'),('*', '*')]


class CDSWindowFeat:
    def __init__(self):
        """
        Input filename: should be a fasta filename containing just CDS sequences (in-frame)
        """
        #self.input_filename = input_filename
        
        cdef string x
        cdef string y
        cdef int k

        self.diamino_range = [1,2,3,4] # parameter k
        self.amino_total = 0
        self.amino_count = {}
        # init amino count
        for x in AMINO_LETTERS:
            self.amino_count[x] = 0

        self.codon_total = 0
        self.codon_count = {}
        for y in CODON_LETTERS:
            self.codon_count[y] = 0

        self.diamino_count = {} # dict of k --> (A_i, B_j) --> count
        self.diamino_total = {}
        for k in self.diamino_range:
            self.diamino_total[k] = 0
            self.diamino_count[k] = {}
            for x, y in DIAMINO_LETTERS:
                self.diamino_count[k][(x, y)] = 0

        self.amino_freq = {}
        self.codon_freq = {}
        self.diamino_freq = {}
        self.diamino_changed = {} # (k, x, y)
#

#
#    def __add__(self, other):
#        assert self.diamino_range == other.diamino_range
#        result = CDSWindowFeat()
#        for x in CDSWindowFeat.AMINO_LETTERS:
#            result.amino_count[x] = self.amino_count[x] + other.amino_count[x]
#        for a in CDSWindowFeat.AMINO_LETTERS:
#            for b in CDSWindowFeat.AMINO_LETTERS:
#                for k in self.diamino_range:
#                    result.diamino_count[k][(a,b)] = self.diamino_count[k][(a,b)] + other.diamino_count[k][(a,b)]
#        for a,b,c in itertools.product(CDSWindowFeat.NT_LETTERS, CDSWindowFeat.NT_LETTERS, CDSWindowFeat.NT_LETTERS):
#            result.codon_count[a+b+c] = self.codon_count[a+b+c] + other.codon_count[a+b+c]
#
#        result.codon_total = self.codon_total + other.codon_total
#        result.amino_total = self.amino_total + other.amino_total
#        for k in self.diamino_range:
#            result.diamino_total[k] = self.diamino_total[k] + other.diamino_total[k]
#        return result


    def calc_amino_count(self, char* aa_seq, int factor=1):
        cdef string x
        for x in aa_seq:
            self.amino_count[x] += 1 * factor
            self.amino_total += 1 * factor

    def calc_codon_count(self, char* nt_seq, int nt_len, int factor=1):
        cdef int i
        for i in xrange(0, nt_len-2, 3):
            self.codon_count[nt_seq[i:i+3]] += 1 * factor
            self.codon_total += 1 * factor

    def calc_diamino_count(self, char* aa_seq, int aa_len):
        cdef int i, k
        cdef char x, y
        for i in xrange(aa_len-1):
            for k in self.diamino_range:
                if i + k < aa_len:
                    x = aa_seq[i]
                    y = aa_seq[i+k]
                    self.diamino_count[k][(chr(x), chr(y))] += 1 
                    self.diamino_total[k] += 1 
                    
    def deduct_diamino_count(self, char* aa_seq, int aa_len, int i_range):
        cdef int i, k
        cdef char x, y
        for i in xrange(i_range):
            for k in self.diamino_range:
                if i + k < aa_len:
                    x = aa_seq[i]
                    y = aa_seq[i+k]
                    self.diamino_count[k][(chr(x), chr(y))] -= 1 
                    self.diamino_total[k] -= 1 
                    self.diamino_changed[(k, chr(x), chr(y))] = 1
                    
    def add_diamino_count(self, char* aa_seq, int aa_len, int i_range):
        cdef int i, k
        cdef char x, y
        for i in xrange(aa_len-1, aa_len-1-i_range, -1):
            for k in self.diamino_range:
                if i - k >= 0:
                    x = aa_seq[i-k]
                    y = aa_seq[i]
                    self.diamino_count[k][(chr(x), chr(y))] += 1 
                    self.diamino_total[k] += 1 
                    self.diamino_changed[(k, chr(x), chr(y))] = 1

    def get_amino_freq(self, object pseudo=None, double alpha=0):
        cdef string x
        cdef double a, b
        self.amino_freq = {}
        for x in AMINO_LETTERS:
            a = self.amino_count[x]
            b = self.amino_total
            if pseudo is not None and alpha > 0:
                a += alpha*pseudo.amino_count[x]
                b += alpha*pseudo.amino_total
            self.amino_freq[x] = a * 1. / b
        return self.amino_freq

    def get_codon_freq(self, object pseudo=None, double alpha=0):
        cdef string x
        cdef double a, b
        self.codon_freq = {}
        for x in CODON_LETTERS:
            a = self.codon_count[x]
            b = self.codon_total
            if pseudo is not None and alpha > 0:
                #print "adding pseudo for codon ferq", alpha
                a += alpha*pseudo.codon_count[x]
                b += alpha*pseudo.codon_total
            self.codon_freq[x] = a * 1. / b
        return self.codon_freq

    def get_diamino_freq(self, object pseudo=None, double alpha=0, bool clear_dict=True):
        cdef int k
        cdef double a, b
        cdef string x, y
        if clear_dict: self.diamino_freq = {}
        for k in self.diamino_range:
            if clear_dict: self.diamino_freq[k] = {}
            for x, y in DIAMINO_LETTERS:
                if not clear_dict and (k, x, y) not in self.diamino_changed: continue
                a = self.diamino_count[k][(x, y)]
                b = self.diamino_total[k]
                if pseudo is not None and alpha> 0:
                    a += alpha*pseudo.diamino_count[k][(x, y)]
                    b += alpha*pseudo.diamino_total[k]
                self.diamino_freq[k][(x,y)] = a * 1. / b
        return self.diamino_freq


#def make_amino_scores(freq):
#    arr = []
#    # make sure it's in the right *order*
#    for x in CDSWindowFeat.AMINO_LETTERS:
#        arr.append(freq[x])
#    return arr
#
#def make_diamino_scores(freq, diamino_range):
#    """
#    In the order of (k, A_i, B_j)
#    """
#    arr = []
#    for k in diamino_range:
#        for a in CDSWindowFeat.AMINO_LETTERS:
#            for b in CDSWindowFeat.AMINO_LETTERS:
#                arr.append(freq[k][(a, b)])
#    return arr
#
#def make_codon_scores(aa_freq, codon_freq):
#    arr = []
#    letters = CDSWindowFeat.NT_LETTERS
#    for codon, aa in CodonTable.standard_dna_table.forward_table.iteritems():
#        arr.append(aa_freq[aa]*1./codon_freq[codon])
#    for codon in CodonTable.standard_dna_table.stop_codons:
#        arr.append(aa_freq[CDSWindowFeat.STOP_AMINO]*1./codon_freq[codon])
#    return arr
#
#import c_ORFscores as sp2
#def make_data(seq, pseudo, window_size=96, step_size=3, frame_shift=0):
#    n = len(seq)
#    data = []
#    for i in xrange(frame_shift, n-window_size, step_size):
#        s = seq[i:i+window_size].tostring()
#        a = seq[i:i+window_size].translate()
#        o = sp2.CDSWindowFeat()
#        o.calc_amino_count(a)
#        o.calc_diamino_count(a)
#        o.calc_codon_count(s)
#        aa_freq = o.get_amino_freq(pseudo, .00001)
#        di_freq = o.get_diamino_freq(pseudo, .00001)
#        codon_freq = o.get_codon_freq(pseudo, .00001)
#        arr = make_amino_scores(aa_freq)+make_diamino_scores(di_freq,o.diamino_range)+make_codon_scores(aa_freq,codon_freq)
#        data.append(arr)
#    return data
#
#def find_start_stop_codons(seq):
#    """
#    Given a sequence, return a dict of:
#      frame (0, 1, 2) --> n-th codon --> "TAA/TAG/TGA"
#      frame (0, 1, 2) --> n-th codon --> "TAG"
#
#    Returns: start_dict, stop_dict
#    """
#    stop_d = {0: {}, 1: {}, 2: {}}
#    start_d = {0: {}, 1: {}, 2: {}}
#    for i in xrange(len(seq)-3):
#        if seq[i:i+3] in CodonTable.standard_dna_table.stop_codons:
#            frame = i % 3
#            nth = (i-frame) / 3
#            stop_d[i%3][nth] = seq[i:i+3]
#        elif seq[i:i+3] in CodonTable.standard_dna_table.start_codons:
#            frame = i % 3
#            nth = (i-frame) / 3
#            start_d[i%3][nth] = seq[i:i+3]
#    return start_d, stop_d
#
#
#def find_indel(seq, max_i, bdt, background, target_frame, old_score):
#    best_score, best_moves = old_score, []
#    # try deleting one base first
#    for i in xrange(max_i):
#        new_seq = seq[:i] + seq[(i+1):]
#        start_dict, stop_dict = find_start_stop_codons(new_seq.tostring())
#        if len(stop_dict[target_frame]) > 0: continue # not feasible if introduces a stop codon
#        stuff = make_data(new_seq, background, frame_shift=target_frame)
#        score = sum(bdt.predict(stuff))
#        if score > best_score:
#            best_score, best_moves = score, [i]
#        elif score == best_score:
#            best_moves.append(i)
#
#    return best_score, best_moves
