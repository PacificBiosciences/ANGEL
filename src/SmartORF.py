__author__ = 'lachesis'

import os, sys, subprocess, time
from cPickle import load, dump
from collections import defaultdict
from multiprocessing import Process, Queue, Pool
import numpy as np
from sklearn.ensemble import AdaBoostClassifier
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Angel
from Angel import c_ORFscores, ORFscores, DumbORF
from Angel.ORFutils import write_CDS_n_PEP, convert_non_ATCG

sys.modules['c_ORFscores'] = Angel.c_ORFscores # need this for unpickling to work

MAX_RECORD_CHUNK = 100




def add_to_background(o_all, records):
    for r in records:
        o_all.calc_codon_count(str(r.seq).upper(), len(r.seq))
        aa = str(r.seq.translate())
        o_all.calc_amino_count(aa, len(aa))
        o_all.calc_diamino_count(aa, len(aa))

def add_data_worker(o_all, records, frames, queue):
    #result = []
    for rec in records:
        for i in frames:
            print >> sys.stderr, "processing record {0}, frame {1}".format(rec.id, i)
            stuff = ORFscores.make_data_smart(rec.seq, o_all, frame_shift=i)
            print >> sys.stderr, "putting into queue", rec.id
            queue.put(stuff)
            #result += stuff
            print >> sys.stderr, "done for ", rec.id
    print >> sys.stderr, "Done with records"
    #with open(filename, 'w') as f:
    #    dump(result, f)

def get_data_parallel(o_all, records, frames, num_workers):
    data = []
    workers = []
    queue = Queue()
    bin_size = len(records) / num_workers + 1
    for i in xrange(num_workers):
        #filename = "tmp.{0}.pickle".format(i)
        #p = Process(target=add_data_worker, args=(o_all, records[bin_size*i:bin_size*(i+1)], frames, queue, filename))
        p = Process(target=add_data_worker, args=(o_all, records[bin_size*i:bin_size*(i+1)], frames, queue))
        workers.append(p)

    print >> sys.stderr, "launching all workers"
    for p in workers:
        print >> sys.stderr, "launching worker", p.name
        p.start()

    time.sleep(10) # simply wait 30 sec
    # print >> sys.stderr, "waiting for workers to finish...."
    total = len(records)*len(frames)
    fail_count = 0
    for i in xrange(total):
        obj = queue.get(timeout=60)
        if obj is not None:
            data += obj
        else:
            print >> sys.stderr, "record {0} of {1} waittime exceeded! give up!".format(i, total)
            fail_count += 1
        if fail_count > 10:
            print >> sys.stderr, "failed waittime more than 10 times. stop!"
            break

    # must purge everything in queue
    while not queue.empty():
        data += queue.get()

    print >> sys.stderr, "purged queue. now can join workers."

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

    # Queue is very inefficient for large data passing
    # instead break the records up into chunk sizes and just combine results together
    print >> sys.stderr, "running get_data_parallel for coding, chunk 0"
    data_pos = get_data_parallel(o_all, coding, [0], num_workers)
    data_neg = get_data_parallel(o_all, utr, [0, 1, 2], num_workers)
#    num_coding = len(coding)
#    data_pos = get_data_parallel(o_all, coding[:MAX_RECORD_CHUNK], [0], num_workers)
#    for i in xrange(1, num_coding/MAX_RECORD_CHUNK + (num_coding%MAX_RECORD_CHUNK>0)):
#        print >> sys.stderr, "running get_data_parallel for coding, chunk", i
#        data_pos += get_data_parallel(o_all, coding[i*MAX_RECORD_CHUNK:(i+1)*MAX_RECORD_CHUNK], [0], num_workers)##
#
#    print >> sys.stderr, "running get_data_parallel for UTR, chunk 0"
#    num_utr = len(utr)
#    data_neg = get_data_parallel(o_all, utr[:MAX_RECORD_CHUNK], [0, 1, 2], num_workers)
#    for i in xrange(1, num_utr/MAX_RECORD_CHUNK + (num_utr%MAX_RECORD_CHUNK>0)):
#        print >> sys.stderr, "running get_data_parallel for UTR, chunk", i
#        data_neg += get_data_parallel(o_all, utr[i*MAX_RECORD_CHUNK:(i+1)*MAX_RECORD_CHUNK], [0, 1, 2], num_workers)

    print >> sys.stderr, "size of neg training data: {0}, pos training data: {1}".format(\
        len(data_neg), len(data_pos))

    print >> sys.stderr, "using first 10,000 training pos/neg only"
    data_neg = data_neg[:10000]
    data_pos = data_pos[:10000]
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


def ANGEL_predict_worker(input_fasta, output_prefix, bdt, o_all, min_ANGEL_aa_length=50, min_dumb_aa_length=100, use_rev_strand=False, output_mode='best', max_angel_secondORF_distance=10, starting_index=1):
    """
    Output Mode is either "best" or "all"

    If "all" and + strand only: ANGEL+, dumb+  (both subject to its length threshold)
    If "all" and - strand also: ANGEL+, dumb+, ANGEL-, dumb-

    If "best" and + strand only: argmax_(ANGEL+, dumb+)
    If "best" and - strand also: argmax_(ANGEL+, dumb+, ANGEL-, dumb-)

    For dumb, pick only the longest ORF ouf of the 3 possible frames for that strand.
    For ANGEL, if there are multiple ORFs (suspicious), the longest one is chosen as the "length" to beat dumb,
              and if ANGEL is chosen as output, all ORFs are output.
    """

    for rec in SeqIO.parse(open(input_fasta), 'fasta'):
        ORFs = []
        # convert any non-ATCG to 'A'
        rec.seq = Seq(convert_non_ATCG(str(rec.seq), replace_with='A'))
        seq_len = len(rec.seq)
        n, m = seq_len/3, seq_len%3
        print >> sys.stderr, "predicting for", rec.id
        # (1a) predict on + strand for ANGEL
        result = defaultdict(lambda: []) # frame --> list of (type, start, end)
        stuff = [] # (frame, type, start, end)  # this should eventually replace result, keeping result for now.
        flag, name, good = ORFscores.predict_ORF(rec, bdt, o_all, min_aa_len=min_ANGEL_aa_length)
        #print >> sys.stderr, flag, name, good
        for _frame, _stop, _start in good:
            s = _start * 3 + _frame if _start is not None else _frame
            e = _stop * 3 + _frame + 3 if _stop is not None else n*3 + (_frame if m >= _frame else 0)
            result[_frame].append((flag, s, e))
            stuff.append((_frame, flag, s, e))

        # REGARDLESS OF FRAME, only keep the first ORF unless the later ones overlap or is sufficiently close
        stuff.sort(key=lambda (a,b,c,d): (c, d-c)) # sort by start, then length
        i = 1
        while i < len(stuff):
            if stuff[i-1][3]-max_angel_secondORF_distance <= stuff[i][2] <= stuff[i-1][3]+max_angel_secondORF_distance:
                i += 1
            else: # is too far, kick it!
                stuff.pop(i)
        # put stuff back into result as a dict
        result = defaultdict(lambda: []) # result is effectively overwritten, in the future I can just remove the result in the lines above
        for _frame, _flag, _start, _end in stuff:
            result[_frame].append((_flag, _start, _end))

        if len(result) > 0:
            ORFs.append((rec, result, '+'))


        # (1b) run dumb ORFs which returns the frame with longest ORF as a dict frame -> (flag,s,e) or None
        dumb = DumbORF.predict_longest_ORFs(str(rec.seq).upper(), min_dumb_aa_length)
        if dumb is not None:
            ORFs.append((rec, dumb, '+'))

        # (2a) see if need to predict on - strand
        #      if need to, create a rec2 that has the rev complement
        if use_rev_strand:
            rec2 = SeqRecord(rec.seq.reverse_complement(), id=rec.id, description=rec.description)
            result = defaultdict(lambda: []) # frame --> list of (type, start, end)
            flag, name, good = ORFscores.predict_ORF(rec2, bdt, o_all, min_aa_len=min_ANGEL_aa_length)
            for _frame, _stop, _start in good:
                s = _start * 3 + _frame if _start is not None else _frame
                e = _stop * 3 + _frame + 3 if _stop is not None else n*3 + (_frame if m >= _frame else 0)
                result[_frame].append((flag, s, e))
            # for each frame, only keep the first ORF unless the later ones overlap or is sufficiently close
            for _frame in result:
                stuff = result[_frame]
                stuff.sort(key=lambda (a,b,c): (b,c-b)) # sort by start, then length
                i = 1
                while i < len(stuff):
                    if stuff[i-1][2]-max_angel_secondORF_distance <= stuff[i][1] <= stuff[i-1][2]+max_angel_secondORF_distance:
                        pass
                    else: # is too far, kick it!
                        stuff.pop(i)
                result[_frame] = stuff

            if len(result) > 0:
                ORFs.append((rec, result, '-')) # NOTE: sending rec instead of rec2 here is CORRECT
            dumb = DumbORF.predict_longest_ORFs(str(rec2.seq).upper(), min_dumb_aa_length)
            if dumb is not None:
                ORFs.append((rec, dumb, '-'))

        # now decide what to output from ORFs
        # if output_mode:all, just output everything
        # if output_mode:best, pick the longest one

        if output_mode == 'best' and len(ORFs)>0:
            #print >> sys.stderr, "output mode: best"
            #print >> sys.stderr, ORFs
            best_rec, best_result, best_strand = ORFs[0]
            best_len = max(max(e-s for (flag,s,e) in v) for v in best_result.itervalues())
            for _rec, _result, _strand in ORFs[1:]:
                _len = max(max(e-s for (flag,s,e) in v) for v in _result.itervalues())
                if _len > best_len:
                    best_rec, best_result, best_strand, best_len = \
                    _rec, _result, _strand, _len
            ORFs = [(best_rec, best_result, best_strand)]
        print >> sys.stderr, "writing result for", rec.id, "to", output_prefix
        #print >> sys.stderr, "current ORFs:", ORFs
        starting_index = write_CDS_n_PEP(ORFs, output_prefix, min_utr_length=50, append_file=True, starting_index=starting_index)
    print >> sys.stderr, "ALL DONE for", output_prefix
    os.system("touch {0}.DONE".format(output_prefix))

def ANGEL_predict_worker_helper(args):
    return ANGEL_predict_worker(*args)

def distribute_ANGEL_predict(fasta_filename, output_prefix, bdt_pickle_filename, num_workers=5, min_ANGEL_aa_length=50, min_dumb_aa_length=100, use_rev_strand=False, output_mode='best', max_angel_secondORF_distance=10):
    tmpdir = "ANGEL.tmp." + str(int(time.time()))
    os.makedirs(tmpdir)

    print >> sys.stderr, "Reading classifer pickle:", bdt_pickle_filename
    with open(bdt_pickle_filename, 'rb') as f:
        a = load(f)
        bdt = a['bdt']
        o_all = a['o_all']

    print >> sys.stderr, "Splitting input into chunks for parallelization...."

    ###
    # split the input into chunks of 1000 (fixed) sequences
    # use a Pool to continously work on each of the input
    ###

    total_seqs = 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'): total_seqs += 1
    num_workers = min(num_workers, total_seqs)  # just in case there are fewer sequences than workers

    if total_seqs < num_workers * 1000:
        num_seqs_per_worker = total_seqs / num_workers
    else:
        num_seqs_per_worker = min(1000, total_seqs) # used to be: total_seqs / num_workers + 1
    num_splits = total_seqs / num_seqs_per_worker + 1
    handles = [open(os.path.join(tmpdir, output_prefix+'.split_'+str(i)+'.fa'), 'w') for i in xrange(num_splits)]
    i = 0
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        f = handles[i/num_seqs_per_worker]
        f.write(">{0}\n{1}\n".format(r.id, r.seq))
        i += 1
    for f in handles: f.close()
    handles = filter(lambda f: os.stat(f.name).st_size > 0, handles)

    list_of_fasta = [ f.name for f in handles ]

    n = num_seqs_per_worker
    #workers = []
    data = []
    for i, input_fasta in enumerate(list_of_fasta):
        print >> sys.stderr, "Pool worker for", input_fasta
        starting_index = i * n + 1

        data.append((input_fasta, input_fasta+'.ANGEL', bdt, o_all, \
                     min_ANGEL_aa_length, min_dumb_aa_length, use_rev_strand,\
                     output_mode, max_angel_secondORF_distance, starting_index))
        #p = Process(target=ANGEL_predict_worker, args=(input_fasta, input_fasta+'.ANGEL', bdt, o_all, min_ANGEL_aa_length, min_dumb_aa_length, use_rev_strand, output_mode, max_angel_secondORF_distance, starting_index))
        #p.start()
        #workers.append(p)

    #for p in workers:
    #    print >> sys.stderr, "waiting for worker", p.name
    #    p.join()

    pool = Pool(num_workers)
    pool.map(ANGEL_predict_worker_helper, data)

    print >> sys.stderr, "Closing Pool...."
    pool.close()
    print >> sys.stderr, "Joining Pool...."
    pool.join()
    print >> sys.stderr, "All workers completed."

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