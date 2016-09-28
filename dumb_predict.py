#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'

from Angel.DumbORF import transdecoder_main

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser("Program that predicting longest ORFs")
    parser.add_argument("fasta_filename", help="Fasta filename to train on")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("--min_aa_length", default=100, type=int, help="Minimum ORF length (default: 100aa)")
    parser.add_argument("--use_firstORF", default=False, action="store_true", help="Use first ORF instead of longest ORF")
    parser.add_argument("--use_rev_strand", default=False, action="store_true", help="Predict on reverse strand as well (default: off)")
    parser.add_argument("--cpus", default=8, type=int, help="Number of CPUs (default: 8)")

    args = parser.parse_args()

    transdecoder_main(args.fasta_filename, args.output_prefix, args.min_aa_length, args.use_rev_strand, args.use_firstORF, args.cpus)
