ANGEL: Robust Open Reading Frame prediction
=====
Last Updated: 6/11/2014


The program is divided into three parts :

-- dumb ORF prediction: ORF prediction that outputs all longest ORFs in all frames. Can be used to create a top training dataset.

-- ANGEL classifier training: Train a coding potential classifier based on given training data. Must provide both CDS and UTR for positive and negative training set. 

-- ANGEL robust ORF prediction: ORF prediction based on both the ANGEL classifier and the dumb ORF prediction.


The ANGEL classifer is based on the classifer described in the Shimizu *et al.*, "**ANGLE: a sequencing errors resistant program for predicting protein coding regions in unfinished cDNA.**", *J Bioinform Comput Biol*, (2006). The naming change to ANGEL is intentional to disambiguate from the author implementation of ANGLE.


## SOFTWARE REQUIREMENT
You need to install [CD-HIT](http://www.bioinformatics.org/cd-hit/) and have it available in your ``$PATH`` variable to run dumb ORF prediction.


## INSTALLATION
We recommend that you set up and activate a virtual environment before installation. See [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Setting-up-virtualenv-and-installing-pbtranscript-tofu) for installation details.

You can download this GitHub repository in many ways. Here we assume you will be using git clone:

```
git clone https://github.com/Magdoll/ANGEL.git
cd ANGEL
python setup.py build
python setup.py install
```


## USAGE

#### Dumb ORF prediction

`dumb_predict.py` takes as input a FASTA file. It outputs all longest ORFs (which could be overlapping) that exceed the user-defined minimum length and has a positive log-odds scores based on hexamer frequencies. 


Usage:

```
dumb_predict.py <fasta_filename> <output_prefix> 
       [--min_aa_length MIN_AA_LENGTH]
       [--use_rev_strand] [--cpus CPUS]
```

See the provided training example below. By default, only the forward strand is used. This is especially true for PacBio transcriptome sequencing output that often already has correct strand.

```
cd ANGEL/training_example
dumb_predict.py test.human_1000seqs.fa test.human.dumb --min_aa_length 300 --cpus 24
```

The output consists of *test.human.dumb.final.pep*, *test.human.dumb.final.cds* and *test.human.dumb.final.utr*, which are the results of longest ORF prediction.

#### Creating a non-redundant training dataset

Redundancy in the training dataset, such as highly identical CDS sequences from alternative isoforms or homologous genes, can skew the classifier training. The script `angel_make_training_set.py` clusters an input set of CDS sequences into non-redundant clusters, and outputs a selective subset for training data.

Note that the training dataset does **not** have to be the same as the input to ANGEL ORF prediction below. You can use any curated dataset like Gencode, RefSeq, or others so long as they are from the same or similar species so that the classifier can be trained properly.

Usage:

```
angel_make_training_set.py <input_prefix> <output_prefix>
     [--use_top USE_TOP] [--random] [--cpus CPUS]
```

Using the output from dumb ORF prediction above, the command is:

```
angel_make_training_set.py test.human.dumb.final test.human.dumb.final.training --random --cpus 24
```

Here we use the `--random` parameter to randomly select non-redundant CDS sequences, instead of choosing the longest CDS sequences.

The output files are: *test.human.dumb.final.training.cds*, *test.human.dumb.final.training.utr* and *test.human.dumb.final.training.pep*.


#### ANGEL classifer training


`angel_train.py` takes a CDS FASTA file and a UTR FASTA file and outputs a trained classifier pickle file. 

Usage:

```
angel_train.py <cds_filename> <utr_filename> <output_pickle> [--cpus CPUS]
```

For example:

```
angel_train.py test.human.dumb.final.training.cds test.human.dumb.final.training.utr \
      test.human.dumb.final.classifier.pickle --cpus 24
```

On a typical 500-sequence training dataset, the training may take several hours. The 500-sequence human MCF-7 training set took 4 hours on a 24-core machine.


#### ANGEL ORF prediction

Usage:

```
angel_predict.py <input_fasta> <classifier_pickle> <output_prefix>
       [--min_angel_aa_length] [--min_dumb_aa_length] 
       [--use_rev_strand] [--output_rev_only_if_longer] 
       [--cpus]
```

For each sequence, ANGEL uses the classifer to predict the coding potential of each codon in each of the three frames, then finds the most likely stretch of window that is the open reading frame. It then does the following:

* If there is only one predicted ORF, tag it as "confident". In this case, it is unlikely there are sequencing errors.
* If there are multiple ORFs but only one exceeds `min_angel_aa_length` threshold, output that one ORF and tag it as "likely". In this case, the program is semi-positive that there are no sequencing errors or that the error occurs at the ends of the CDS, allowing successful prediction of a relatively long continuous ORF.
* If there are multiple ORFs (could be in multiple frames) exceeding the length threshold, output all of them and tag them as "suspicious". In this case, the ORFs could be genuine complete ORFs (indicating a polycistronic transcript or alternative ORFs) or fragments of the same ORF that has frame shift due to uncorrected sequencing errors.

The longest ORF length from the ANGEL process is recorded as *T*. Then, the same dumb ORF prediction (which simply looks for the longest stretch of ORF without stop codon interruption) is done on the forward three frames. Each predicted dumb ORF is also outputted if, and only if, its length is greater than both *T* and `min_dumb_aa_length`. This is a fallback process in case the ANGEL classifier failed to detect coding potential in the CDS region.

If `--use_rev_strand` is used, then the same process is repeated on the reverse-complement of the sequence. 
If `--output_rev_only_if_longer` is used, then the reverse strand ORF is output only if it is longer than the longest ORF from the forward strand.


For example:

```
angel_predict.py test.human_1000seqs.fa human.MCF7.random_500_for_training.pickle test.human \
      --use_rev_strand --output_rev_only_if_longer --min_dumb_aa_length 100
```

The output files are: <output_prefix>.ANGEL.cds, <output_prefix>.ANGEL.pep, <output_prefix>.ANGEL.utr.

Each output sequence ID has the format:
```
<seq_id> type:<tag>-<completeness> len:<ORF length (aa)> strand:<strand> pos:<CDS range>
```

Where *tag* is "confident", "likely", or "suspicious" for ANGEL predictions, and "dumb" for dumb ORF predictions.
*completeness* is either "complete", "5partial", "3partial", or "internal" based on the presence or absence of start and stop codons.


## LICENSE
```
#################################################################################$$
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$
```
