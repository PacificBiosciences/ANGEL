#cython: language_level=3
__author__ = 'etseng@pacificbiosciences.com'

import itertools
from Bio.Alphabet import IUPAC
import Bio.Data.CodonTable as CodonTable

from libcpp.utility cimport pair
from libcpp.map cimport map 
from libcpp.string cimport string
from libcpp cimport bool

AMINO_LETTERS = [b'A', b'C', b'D', b'E', b'F', b'G',\
		b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R',\
		b'S', b'T', b'V', b'W', b'Y', b'*']
NT_LETTERS = [b'G', b'A', b'T', b'C']
CODON_LETTERS = [b"GGG",b"GGA",b"GGT",b"GGC",b"GAG",b"GAA",b"GAT",b"GAC",b"GTG",b"GTA",b"GTT",b"GTC",b"GCG",b"GCA",b"GCT",b"GCC",b"AGG",b"AGA",b"AGT",b"AGC",b"AAG",b"AAA",b"AAT",b"AAC",b"ATG",b"ATA",b"ATT",b"ATC",b"ACG",b"ACA",b"ACT",b"ACC",b"TGG",b"TGA",b"TGT",b"TGC",b"TAG",b"TAA",b"TAT",b"TAC",b"TTG",b"TTA",b"TTT",b"TTC",b"TCG",b"TCA",b"TCT",b"TCC",b"CGG",b"CGA",b"CGT",b"CGC",b"CAG",b"CAA",b"CAT",b"CAC",b"CTG",b"CTA",b"CTT",b"CTC",b"CCG",b"CCA",b"CCT",b"CCC"]
DIAMINO_LETTERS = [\
    (b'A', b'A'),(b'A', b'C'),(b'A', b'D'),(b'A', b'E'),(b'A', b'F'),(b'A', b'G'),(b'A', b'H'),(b'A', b'I'),(b'A', b'K'),(b'A', b'L'),(b'A', b'M'),(b'A', b'N'),(b'A', b'P'),(b'A', b'Q'),(b'A', b'R'),(b'A', b'S'),(b'A', b'T'),(b'A', b'V'),(b'A', b'W'),(b'A', b'Y'),(b'A', b'*'),(b'C', b'A'),(b'C', b'C'),(b'C', b'D'),(b'C', b'E'),(b'C', b'F'),(b'C', b'G'),(b'C', b'H'),(b'C', b'I'),(b'C', b'K'),(b'C', b'L'),(b'C', b'M'),(b'C', b'N'),(b'C', b'P'),(b'C', b'Q'),(b'C', b'R'),(b'C', b'S'),(b'C', b'T'),(b'C', b'V'),(b'C', b'W'),(b'C', b'Y'),(b'C', b'*'),(b'D', b'A'),(b'D', b'C'),(b'D', b'D'),(b'D', b'E'),(b'D', b'F'),(b'D', b'G'),(b'D', b'H'),(b'D', b'I'),(b'D', b'K'),(b'D', b'L'),(b'D', b'M'),(b'D', b'N'),(b'D', b'P'),(b'D', b'Q'),(b'D', b'R'),(b'D', b'S'),(b'D', b'T'),(b'D', b'V'),(b'D', b'W'),(b'D', b'Y'),(b'D', b'*'),(b'E', b'A'),(b'E', b'C'),(b'E', b'D'),(b'E', b'E'),(b'E', b'F'),(b'E', b'G'),(b'E', b'H'),(b'E', b'I'),(b'E', b'K'),(b'E', b'L'),(b'E', b'M'),(b'E', b'N'),(b'E', b'P'),(b'E', b'Q'),(b'E', b'R'),(b'E', b'S'),(b'E', b'T'),(b'E', b'V'),(b'E', b'W'),(b'E', b'Y'),(b'E', b'*'),(b'F', b'A'),(b'F', b'C'),(b'F', b'D'),(b'F', b'E'),(b'F', b'F'),(b'F', b'G'),(b'F', b'H'),\
    (b'F', b'I'),(b'F', b'K'),(b'F', b'L'),(b'F', b'M'),(b'F', b'N'),(b'F', b'P'),(b'F', b'Q'),(b'F', b'R'),(b'F', b'S'),(b'F', b'T'),(b'F', b'V'),(b'F', b'W'),(b'F', b'Y'),(b'F', b'*'),(b'G', b'A'),(b'G', b'C'),(b'G', b'D'),(b'G', b'E'),(b'G', b'F'),(b'G', b'G'),(b'G', b'H'),(b'G', b'I'),(b'G', b'K'),(b'G', b'L'),(b'G', b'M'),(b'G', b'N'),(b'G', b'P'),(b'G', b'Q'),(b'G', b'R'),(b'G', b'S'),(b'G', b'T'),(b'G', b'V'),(b'G', b'W'),(b'G', b'Y'),(b'G', b'*'),(b'H', b'A'),(b'H', b'C'),(b'H', b'D'),(b'H', b'E'),(b'H', b'F'),(b'H', b'G'),(b'H', b'H'),(b'H', b'I'),(b'H', b'K'),(b'H', b'L'),(b'H', b'M'),(b'H', b'N'),(b'H', b'P'),(b'H', b'Q'),(b'H', b'R'),(b'H', b'S'),(b'H', b'T'),(b'H', b'V'),(b'H', b'W'),(b'H', b'Y'),(b'H', b'*'),(b'I', b'A'),(b'I', b'C'),(b'I', b'D'),(b'I', b'E'),(b'I', b'F'),(b'I', b'G'),(b'I', b'H'),(b'I', b'I'),(b'I', b'K'),(b'I', b'L'),(b'I', b'M'),(b'I', b'N'),(b'I', b'P'),(b'I', b'Q'),(b'I', b'R'),(b'I', b'S'),(b'I', b'T'),(b'I', b'V'),(b'I', b'W'),(b'I', b'Y'),(b'I', b'*'),(b'K', b'A'),(b'K', b'C'),(b'K', b'D'),(b'K', b'E'),(b'K', b'F'),(b'K', b'G'),(b'K', b'H'),(b'K', b'I'),(b'K', b'K'),(b'K', b'L'),(b'K', b'M'),(b'K', b'N'),(b'K', b'P'),(b'K', b'Q'),\
    (b'K', b'R'),(b'K', b'S'),(b'K', b'T'),(b'K', b'V'),(b'K', b'W'),(b'K', b'Y'),(b'K', b'*'),(b'L', b'A'),(b'L', b'C'),(b'L', b'D'),(b'L', b'E'),(b'L', b'F'),(b'L', b'G'),(b'L', b'H'),(b'L', b'I'),(b'L', b'K'),(b'L', b'L'),(b'L', b'M'),(b'L', b'N'),(b'L', b'P'),(b'L', b'Q'),(b'L', b'R'),(b'L', b'S'),(b'L', b'T'),(b'L', b'V'),(b'L', b'W'),(b'L', b'Y'),(b'L', b'*'),(b'M', b'A'),(b'M', b'C'),(b'M', b'D'),(b'M', b'E'),(b'M', b'F'),(b'M', b'G'),(b'M', b'H'),(b'M', b'I'),(b'M', b'K'),(b'M', b'L'),(b'M', b'M'),(b'M', b'N'),(b'M', b'P'),(b'M', b'Q'),(b'M', b'R'),(b'M', b'S'),(b'M', b'T'),(b'M', b'V'),(b'M', b'W'),(b'M', b'Y'),(b'M', b'*'),(b'N', b'A'),(b'N', b'C'),(b'N', b'D'),(b'N', b'E'),(b'N', b'F'),(b'N', b'G'),(b'N', b'H'),(b'N', b'I'),(b'N', b'K'),(b'N', b'L'),(b'N', b'M'),(b'N', b'N'),(b'N', b'P'),(b'N', b'Q'),(b'N', b'R'),(b'N', b'S'),(b'N', b'T'),(b'N', b'V'),(b'N', b'W'),(b'N', b'Y'),(b'N', b'*'),(b'P', b'A'),(b'P', b'C'),(b'P', b'D'),(b'P', b'E'),(b'P', b'F'),(b'P', b'G'),(b'P', b'H'),(b'P', b'I'),(b'P', b'K'),(b'P', b'L'),(b'P', b'M'),(b'P', b'N'),(b'P', b'P'),(b'P', b'Q'),(b'P', b'R'),(b'P', b'S'),(b'P', b'T'),(b'P', b'V'),(b'P', b'W'),(b'P', b'Y'),(b'P', b'*'),(b'Q', b'A'),(b'Q', b'C'),(b'Q', b'D'),(b'Q', b'E'),(b'Q', b'F'),(b'Q', b'G'),(b'Q', b'H'),(b'Q', b'I'),(b'Q', b'K'),(b'Q', b'L'),(b'Q', b'M'),(b'Q', b'N'),(b'Q', b'P'),(b'Q', b'Q'),(b'Q', b'R'),(b'Q', b'S'),(b'Q', b'T'),(b'Q', b'V'),(b'Q', b'W'),(b'Q', b'Y'),(b'Q', b'*'),(b'R', b'A'),(b'R', b'C'),(b'R', b'D'),(b'R', b'E'),(b'R', b'F'),(b'R', b'G'),(b'R', b'H'),(b'R', b'I'),(b'R', b'K'),(b'R', b'L'),(b'R', b'M'),(b'R', b'N'),(b'R', b'P'),(b'R', b'Q'),(b'R', b'R'),(b'R', b'S'),(b'R', b'T'),(b'R', b'V'),(b'R', b'W'),(b'R', b'Y'),(b'R', b'*'),(b'S', b'A'),(b'S', b'C'),(b'S', b'D'),(b'S', b'E'),(b'S', b'F'),(b'S', b'G'),(b'S', b'H'),(b'S', b'I'),(b'S', b'K'),(b'S', b'L'),(b'S', b'M'),(b'S', b'N'),(b'S', b'P'),(b'S', b'Q'),(b'S', b'R'),(b'S', b'S'),(b'S', b'T'),(b'S', b'V'),(b'S', b'W'),(b'S', b'Y'),(b'S', b'*'),(b'T', b'A'),(b'T', b'C'),(b'T', b'D'),(b'T', b'E'),(b'T', b'F'),(b'T', b'G'),(b'T', b'H'),(b'T', b'I'),(b'T', b'K'),(b'T', b'L'),(b'T', b'M'),(b'T', b'N'),(b'T', b'P'),(b'T', b'Q'),(b'T', b'R'),(b'T', b'S'),(b'T', b'T'),(b'T', b'V'),(b'T', b'W'),(b'T', b'Y'),(b'T', b'*'),(b'V', b'A'),(b'V', b'C'),(b'V', b'D'),(b'V', b'E'),(b'V', b'F'),(b'V', b'G'),(b'V', b'H'),(b'V', b'I'),(b'V', b'K'),(b'V', b'L'),(b'V', b'M'),(b'V', b'N'),(b'V', b'P'),(b'V', b'Q'),(b'V', b'R'),(b'V', b'S'),(b'V', b'T'),(b'V', b'V'),(b'V', b'W'),(b'V', b'Y'),(b'V', b'*'),(b'W', b'A'),(b'W', b'C'),(b'W', b'D'),(b'W', b'E'),(b'W', b'F'),(b'W', b'G'),(b'W', b'H'),(b'W', b'I'),(b'W', b'K'),(b'W', b'L'),(b'W', b'M'),(b'W', b'N'),(b'W', b'P'),(b'W', b'Q'),(b'W', b'R'),(b'W', b'S'),(b'W', b'T'),(b'W', b'V'),(b'W', b'W'),(b'W', b'Y'),(b'W', b'*'),(b'Y', b'A'),(b'Y', b'C'),(b'Y', b'D'),(b'Y', b'E'),(b'Y', b'F'),(b'Y', b'G'),(b'Y', b'H'),(b'Y', b'I'),(b'Y', b'K'),(b'Y', b'L'),(b'Y', b'M'),(b'Y', b'N'),(b'Y', b'P'),(b'Y', b'Q'),(b'Y', b'R'),(b'Y', b'S'),(b'Y', b'T'),(b'Y', b'V'),(b'Y', b'W'),(b'Y', b'Y'),(b'Y', b'*'),(b'*', b'A'),(b'*', b'C'),(b'*', b'D'),(b'*', b'E'),(b'*', b'F'),(b'*', b'G'),(b'*', b'H'),(b'*', b'I'),(b'*', b'K'),(b'*', b'L'),(b'*', b'M'),(b'*', b'N'),(b'*', b'P'),(b'*', b'Q'),(b'*', b'R'),(b'*', b'S'),(b'*', b'T'),(b'*', b'V'),(b'*', b'W'),(b'*', b'Y'),(b'*', b'*')]

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
        cdef bytes x
        cdef bytes aa_seq2 = aa_seq
        for x in aa_seq2:
            self.amino_count[x] += 1 * factor
            self.amino_total += 1 * factor

    def calc_codon_count(self, char* nt_seq, int nt_len, int factor=1):
        cdef int i
        for i in xrange(0, nt_len-2, 3):
            self.codon_count[nt_seq[i:i+3]] += 1 * factor
            self.codon_total += 1 * factor

    def calc_diamino_count(self, char* aa_seq, int aa_len):
        cdef int i, k
        cdef bytes x, y
        for i in xrange(aa_len-1):
            for k in self.diamino_range:
                if i + k < aa_len:
                    x = aa_seq[i]
                    y = aa_seq[i+k]
                    self.diamino_count[k][(x,y)] += 1 
                    self.diamino_total[k] += 1 
                    
    def deduct_diamino_count(self, char* aa_seq, int aa_len, int i_range):
        cdef int i, k
        cdef bytes x, y
        for i in xrange(i_range):
            for k in self.diamino_range:
                if i + k < aa_len:
                    x = aa_seq[i]
                    y = aa_seq[i+k]
                    self.diamino_count[k][(x, y)] -= 1 
                    self.diamino_total[k] -= 1 
                    self.diamino_changed[(k, x, y)] = 1
                    
    def add_diamino_count(self, char* aa_seq, int aa_len, int i_range):
        cdef int i, k
        cdef bytes x, y
        for i in xrange(aa_len-1, aa_len-1-i_range, -1):
            for k in self.diamino_range:
                if i - k >= 0:
                    x = aa_seq[i-k]
                    y = aa_seq[i]
                    self.diamino_count[k][(x, y)] += 1 
                    self.diamino_total[k] += 1 
                    self.diamino_changed[(k, x, y)] = 1

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
