ANGEL: Robust Open Reading Frame prediction
=====

The program is divided into three parts:

-- dumb ORF prediction: ORF prediction that outputs all longest ORFs in all frames. Can be used to create a top training dataset.

-- ANGEL classifier training: train a coding potential classifier based on given training data. Must provide both CDS and UTR for positive and negative training set. 

-- ANGEL robust ORF prediction: ORF prediction based on both the ANGEL classifier and the dumb ORF prediction.


The ANGEL classifer is based on the classifer described in the Shimizu *et al.*, "**ANGLE: a sequencing errors resistant program for predicting protein coding regions in unfinished cDNA.**", *J Bioinform Comput Biol*, (2006). The naming change to ANGEL is intentional to disambiguate from the author implementation of ANGLE.


## SOFTWARE REQUIREMENT
You will need to install [CD-HIT](http://www.bioinformatics.org/cd-hit/) and have it available in your $PATH variable to run dumb ORF prediction.


## INSTALLATION
It is recommended that you set up and activate a virtual environment before installation. See [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Setting-up-virtualenv-and-installing-pbtranscript-tofu) howto.

You can download this GitHub repository in many ways, here we assume you will be using git clone:

```
git clone https://github.com/Magdoll/ANGEL.git
cd ANGEL
python setup.py build
python setup.py install
```


## USAGE

#### Dumb ORF prediction

`dumb_predict.py` takes as input a fasta file and outputs all longest ORFs (could be overlapping) that exceed the user-defined minimum length. In addition, it uses [CD-HIT](http://www.bioinformatics.org/cd-hit/) to remove redundancy to create a "top training set" for ANGEL classifier training. 


The usage is:

```
dumb_predict.py fasta_filename output_prefix 
       [--use_top USE_TOP] [--min_aa_length MIN_AA_LENGTH]
       [--use_rev_strand] [--cpus CPUS]
```

We will use the provided training example below. For runtime purpose, we set `--use_top` to 50 here but for proper training we would recommend setting to 500. By default, only the forward strand is used. This is especially true for PacBio transcriptome sequencing output that often already has correct strand.

```
cd ANGEL/training_example
dumb_predict.py test.human_1000seqs.fa test.human --use_top 50 --min_aa_length 300 --cpus 24
```

The output consists of:

test.human.pep, test.human.cds, test.human.utr: which are the results of longest ORF prediction
test.human.training_50.cds, test.human.training_50.utr: which are the top 50 non-redundant CDS/UTR sequences that can be used as training data for ANGEL classifier 


#### ANGEL classifer training

`angel_train.py` takes a CDS fasta file and a UTR fasta file and outputs a trained classifier pickle file. Note that the training dataset does not have to be the same as the input to ANGEL ORF prediction below. You can use any curated dataset like Gencode, RefSeq, or another dataset so long as they are from the same or similar species so the classifier can be trained properly.

The usage is:

```
angel_train.py <cds_filename> <utr_filename> <output_pickle> [--cpus CPUS]
```

For example:

```
cd ANGEL/training_example
angel_train.py test.human.training_50.cds test.human.training_50.utr test.human.training_50.classifier.pickle --cpus 24
```

#### ANGEL ORF prediction



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
