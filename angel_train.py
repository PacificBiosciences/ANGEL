#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'

from Angel.SmartORF import ANGEL_training

if __name__ == "__main__":
    #ANGEL_training(cds_filename, utr_filename, output_pickle, num_workers=3
    from argparse import ArgumentParser
    parser = ArgumentParser("Program for training ANGEL classifier")
    parser.add_argument("cds_filename", help="CDS fasta filename to train on")
    parser.add_argument("utr_filename", help="UTR fasta filename to train on")
    parser.add_argument("output_pickle", help="Output pickle filename")
    parser.add_argument("--cpus", default=8, type=int, help="Number of CPUs (default: 8)")

    args = parser.parse_args()

    ANGEL_training(args.cds_filename, args.utr_filename, args.output_pickle, num_workers=args.cpus)