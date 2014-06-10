__author__ = 'lachesis'

import os, sys, subprocess, time
from cPickle import load, dump
from collections import defaultdict
from multiprocessing import Process, Queue
import numpy as np
from sklearn.ensemble import AdaBoostClassifier
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Angel
from Angel import c_ORFscores, ORFscores, DumbORF
from Angel.ORFutils import write_CDS_n_PEP

sys.modules['c_ORFscores'] = Angel.c_ORFscores # need this for unpickling to work

def add_to_background(o_all, records):
    for r in records:
        o_all.calc_codon_count(r.seq.tostring().upper(), len(r.seq))
        aa = r.seq.translate().tostring()
        o_all.calc_amino_count(aa, len(aa))
        o_all.calc_diamino_count(aa, len(aa))

def add_data_worker(o_all, records, frames, queue):
    for rec in records:
        for i in frames:
            stuff = ORFscores.make_data_smart(rec.seq, o_all, frame_shift=i)
            print >> sys.stderr, "putting into queue", rec.id
            queue.put(stuff)
            print >> sys.stderr, "done for ", rec.id
    print >> sys.stderr, "Done with records"

def get_data_parallel(o_all, records, frames, num_workers):
    data = []
    workers = []
    queue = Queue()
    bin_size = len(records) / num_workers + 1
    for i in xrange(num_workers):
        p = Process(target=add_data_worker, args=(o_all, records[bin_size*i:bin_size*(i+1)], frames, queue))
        p.start()
        workers.append(p)

    print >> sys.stderr, "waiting for workers to finish...."
    for i in xrange(len(records)*len(frames)):
        data += queue.get()
    for p in workers:
        p.join()
    print >> sys.stderr, "all workers done!"
    print >> sys.stderr, "retrieved {0} data items".format(len(data))
    return data


def ANGEL_training(cds_filename, utr_filename, output_pickle, num_workers=3):
    coding = [ r for r in SeqIO.parse(open(cds_filename), 'fasta') ]
    utr = [ r for r in SeqIO.parse(open(utr_filename), 'fasta') ]

    o_all = c_ORFscores.CDSWindowFeat()
    add_to_background(o_all, coding)
    add_to_background(o_all, utr)

    data_pos = get_data_parallel(o_all, coding, [0], num_workers)
    data_neg = get_data_parallel(o_all, utr, [0, 1, 2], num_workers)

    data = data_neg + data_pos
    target = [0]*len(data_neg) + [1]*len(data_pos)
    data = np.array(data)

    print >> sys.stderr, "data prep done, running classifier...."
    bdt = AdaBoostClassifier(n_estimators=50)
    bdt.fit(data, target)

    print >> sys.stderr, "classifier trained. putting pickle to", output_pickle

    with open(output_pickle, 'wb') as f:
        dump({'bdt':bdt, 'o_all':o_all}, f)

    return data, target, bdt


def ANGEL_predict_worker(input_fasta, output_prefix, bdt, o_all, min_ANGEL_aa_length=50, min_dumb_aa_length=100, use_rev_strand=False, output_rev_only_if_longer=False, starting_index=1):
    for rec in SeqIO.parse(open(input_fasta), 'fasta'):
        ORFs = []
        seq_len = len(rec.seq)
        n, m = len(rec.seq)/3, len(rec.seq)%3
        print >> sys.stderr, "predicting for", rec.id
        # (1a) predict on + strand
        result = defaultdict(lambda: []) # frame --> list of (type, start, end)
        max_angle_predicted_orf_len = min_dumb_aa_length
        flag, name, good = ORFscores.predict_ORF(rec, bdt, o_all, min_aa_len=min_ANGEL_aa_length)
        #print >> sys.stderr, flag, name, good
        for _frame, _stop, _start in good:
            s = _start * 3 + _frame if _start is not None else _frame
            e = _stop * 3 + _frame + 3 if _stop is not None else n*3 + (_frame if m >= _frame else 0)
            result[_frame].append((flag, s, e))
            max_angle_predicted_orf_len = max(max_angle_predicted_orf_len, (e - s)/3 + 1)
        ORFs.append((rec, result, '+'))
        # (1b) run dumb ORFs, if better than longest of ANGEL's output it as well
        dumb = DumbORF.predict_longest_ORFs(rec.seq.tostring().upper(), max_angle_predicted_orf_len)

        if sum(len(v) for v in dumb.itervalues()) > 0:
            ORFs.append((rec, dumb, '+'))
            for v in dumb.itervalues():
                if len(v) > 0:
                    for _flag, _s, _e in v:
                        max_angle_predicted_orf_len = max(max_angle_predicted_orf_len, (_e - _s)/3 + 1)

        # (2a) see if need to predict on - strand
        #      if need to, create a rec2 that has the rev complement
        if use_rev_strand:
            #print "output_rev_only_if_longer:", output_rev_only_if_longer
            if output_rev_only_if_longer: # min aa length must be longer than the forward strand longest prediction
                min_dumb_aa_length_for_rev = max_angle_predicted_orf_len
                min_ANGEL_aa_length_for_rev = max(max_angle_predicted_orf_len, min_ANGEL_aa_length)
            else:
                min_dumb_aa_length_for_rev = min_dumb_aa_length
                min_ANGEL_aa_length_for_rev = min_ANGEL_aa_length
                #print min_dumb_aa_length, min_ANGEL_aa_length
            rec2 = SeqRecord(rec.seq.reverse_complement(), id=rec.id, description=rec.description)
            result = defaultdict(lambda: []) # frame --> list of (type, start, end)
            max_angle_predicted_orf_len = min_dumb_aa_length_for_rev
            #print "calling rev with min_aa_len", min_ANGEL_aa_length
            flag, name, good = ORFscores.predict_ORF(rec2, bdt, o_all, min_aa_len=min_ANGEL_aa_length_for_rev)
            for _frame, _stop, _start in good:
                s = _start * 3 + _frame if _start is not None else _frame
                e = _stop * 3 + _frame + 3 if _stop is not None else n*3 + (_frame if m >= _frame else 0)
                result[_frame].append((flag, s, e))
                max_angle_predicted_orf_len = max(max_angle_predicted_orf_len, (e-s)/3+1)
            ORFs.append((rec, result, '-')) # NOTE: sending rec instead of rec2 here is CORRECT
            dumb = DumbORF.predict_longest_ORFs(rec2.seq.tostring().upper(), max_angle_predicted_orf_len)
            if sum(len(v) for v in dumb.itervalues()) > 0:
                ORFs.append((rec, dumb, '-')) # NOTE: sending rec instead of rec2 here is CORRECT

        starting_index = write_CDS_n_PEP(ORFs, output_prefix, min_utr_length=50, append_file=True, starting_index=starting_index)


def distribute_ANGEL_predict(fasta_filename, output_prefix, bdt_pickle_filename, num_workers=5, min_ANGEL_aa_length=50, min_dumb_aa_length=100, use_rev_strand=False, output_rev_only_if_longer=False):
    tmpdir = "ANGEL.tmp." + str(int(time.time()))
    os.makedirs(tmpdir)

    print >> sys.stderr, "Reading classifer pickle:", bdt_pickle_filename
    with open(bdt_pickle_filename, 'rb') as f:
        a = load(f)
        bdt = a['bdt']
        o_all = a['o_all']

    print >> sys.stderr, "Splitting input into chunks for parallelization...."
    total_seqs = 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'): total_seqs += 1
    num_workers = min(num_workers, total_seqs)
    num_seqs_per_worker = total_seqs / num_workers + 1
    handles = [open(os.path.join(tmpdir, output_prefix+'.split_'+str(i)+'.fa'), 'w') for i in xrange(num_workers)]
    i = 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        f = handles[i/num_seqs_per_worker]
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        i += 1
    for f in handles: f.close()
    handles = filter(lambda f: os.stat(f.name).st_size > 0, handles)

    list_of_fasta = [ f.name for f in handles ]

    n = ((i / num_workers + 1) / 10) * 10
    workers = []
    for i, input_fasta in enumerate(list_of_fasta):
        print >> sys.stderr, "Pool worker for", input_fasta
        starting_index = i * n + 1
        p = Process(target=ANGEL_predict_worker, args=(input_fasta, input_fasta+'.ANGEL', bdt, o_all, min_ANGEL_aa_length, min_dumb_aa_length, use_rev_strand, output_rev_only_if_longer, starting_index))
        p.start()
        workers.append(p)

    for p in workers:
        print >> sys.stderr, "waiting for worker", p.name
        p.join()

    cmd = "cat {0} > {1}.ANGEL.cds".format(" ".join(x+'.ANGEL.cds' for x in list_of_fasta), output_prefix)
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "Trouble running command", cmd
        sys.exit(-1)
    cmd = "cat {0} > {1}.ANGEL.pep".format(" ".join(x+'.ANGEL.pep' for x in list_of_fasta), output_prefix)
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "Trouble running command", cmd
        sys.exit(-1)
    cmd = "cat {0} > {1}.ANGEL.utr".format(" ".join(x+'.ANGEL.utr' for x in list_of_fasta), output_prefix)
    if subprocess.check_call(cmd, shell=True)!=0:
        print >> sys.stderr, "Trouble running command", cmd
        sys.exit(-1)

    print >> sys.stderr, "Output written to {0}.ANGEL.cds, {0}.ANGEL.pep, {0}.ANGEL.utr".format(output_prefix)

#    for x in list_of_fasta:
#        os.remove(x)
#        os.remove(x + '.ANGEL.cds')
#        os.remove(x + '.ANGEL.pep')
#        os.remove(x + '.ANGEL.utr')

#    os.removedirs(tmpdir)