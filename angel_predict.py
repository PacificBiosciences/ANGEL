#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'

from Angel.SmartORF import distribute_ANGEL_predict


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Program for predicting ORFs using trained classifier")
    parser.add_argument("fasta_filename", help="Fasta filename to train on")
    parser.add_argument("classifier_pickle", help="Trained classifier pickle name")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("--min_angel_aa_length", default=50, type=int, help="Minimum ORF length predicted from ANGEL (default: 50 aa)")
    parser.add_argument("--max_angel_secondORF_distance", default=10, type=int, help="Maximum distance allowed for the second ORF from the first ORF to be reported (default: 10 bp)")
    parser.add_argument("--min_dumb_aa_length", default=100, type=int, help="Minimum ORF length predicted from dumbORF (default: 100 aa)")
    parser.add_argument("--use_rev_strand", default=False, action="store_true", help="Predict on reverse strand as well (default: off)")
    parser.add_argument("--output_mode", default="best", choices=["best", "all"], help="Output mode (default: best)")
    parser.add_argument("--cpus", default=8, type=int, help="Number of CPUs (default: 8)")

    args = parser.parse_args()

    distribute_ANGEL_predict(args.fasta_filename, args.output_prefix, args.classifier_pickle, args.cpus, args.min_angel_aa_length, args.min_dumb_aa_length, args.use_rev_strand, args.output_mode, args.max_angel_secondORF_distance)


