#!/usr/bin/env python
__author__ = 'etseng@pacificbiosciences.com'

from Angel.DumbORF import select_for_training

if __name__ == "__main__":
    #select_for_training(input_prefix, output_prefix, use_top=500, random=False, cpus=8)
    from argparse import ArgumentParser
    parser = ArgumentParser("Make training dataset")
    parser.add_argument("input_prefix", help="Input prefix (must have .cds, .utr, .pep)")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument("--use_top", default=500, type=int, help="Size of output training set (default: 500)")
    parser.add_argument("--random", default=False, action="store_true", help="Random selection instead of picking longest ORFs")
    parser.add_argument("--cpus", default=8, type=int, help="Number of CPUs (default: 8)")

    args = parser.parse_args()

    select_for_training(args.input_prefix, args.output_prefix, args.use_top, args.random, args.cpus)