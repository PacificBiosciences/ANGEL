ANGEL: Robust Open Reading Frame prediction
=====

**Latest News**: 
(05.09.2016) Updated instructions on installing dependencies. Now recommend using Anaconda.

(06.27.2014)  See this talk on validating PacBio transcript-based ORF predictions with mass spec data! ANGEL was used for creating the set of ORF predictions from the MCF-7 dataset. [G. Sheynkman: Building a "perfect" proteomics database using PacBio MCF-7 transcriptome data](https://vimeo.com/99358676)


Last Updated: 05/23/2017

Current version: 2.4

===

05.23.2017 New in version 2.4

* An user pointed out a bug in SmartORF.py:169 that would cause an infinite loop if `--use_rev_strand` is called in angel_train.py. Fixed!


11.08.2016 New in version 2.3

* MINOR bug in choosing "best" ORF across frames fixed. v2.2 users recommended update to v2.3 for cleaner results.


11.03.2016 New in version 2.2

* FATAL bug in inplementing max_angel_secondORF_distance fixed!!! v2.1 users must update to v2.2.

09.28.2016 New in version 2.1
* added `--use_firstORF` option in `dumb_predict.py` that outputs first ORF instead of longest ORF
* added `--max_angel_secondORF_distance` option in `angel_predict.py` that only outputs later ORFs if they are close enough to the previous ORFs

===

* <a href="#install">Installation

* <a href="#dumb">Dumb ORF prediction</a>
* <a href="#trainset">Creating a non-redundant training dataset</a>
* <a href="#training">ANGEL classifier training</a>
* <a href="#angelpredict">ANGEL ORF prediction</a>
* <a href="#errgenome">Using genomic to error correct first before ORF prediction</a>

===

The program is divided into three parts :

* **Dumb ORF prediction**: ORF prediction that outputs the longest ORF (or first ORF) in all frames. Can be used to create a top training dataset.

* **ANGEL classifier training**: Train a coding potential classifier based on given training data. Must provide both CDS and UTR for positive and negative training set. 

* **ANGEL robust ORF prediction**: ORF prediction based on both the ANGEL classifier and the dumb ORF prediction.


The ANGEL classifier is based on the classifier described in Shimizu *et al.*, "**ANGLE: a sequencing errors resistant program for predicting protein coding regions in unfinished cDNA.**", *J Bioinform Comput Biol*, (2006). The naming change to ANGEL is intentional to differentiate it from the author implementation of ANGLE.

**NOTE**: for both dumb and ANGEL ORF prediction, it is recommended that the input sequences have at least 99% accuracy. This means either short read assembled transcripts or for PacBio, the output from running the [Iso-Seq pipeline](https://github.com/PacificBiosciences/cDNA_primer/) which are Quiver-polished high-quality consensus sequences. ANGEL has not been tested on PacBio subread-level or ReadsOfInsert-level sequences.


<a href="install"/>
## SOFTWARE REQUIREMENT
You need to install [CD-HIT](http://www.bioinformatics.org/cd-hit/) and have it available in your ``$PATH`` variable to run dumb ORF prediction.

## PYTHON PREREQUISITE

The python dependencies for ANGEL are:
* numpy
* Biopython
* scikit-learn


You can install them using several ways, but the most recommended one is using package manager Anaconda (similar to what is done for [Cogent](https://github.com/Magdoll/Cogent/wiki/Installing-Cogent)).

In fact, if you already have [Cogent](https://github.com/Magdoll/Cogent/wiki/Installing-Cogent) installed, the only thing you need to do is activate the environment and install `scikit-learn`.

### Option 1: Install ANGEL using Anaconda

We will first install Anaconda and create a virtualenv under it. The install directions are the same as those for Cogent, so I'm naming the environment `AnaCogent` as well.


(1) [Download](https://www.continuum.io/downloads) Anaconda latest version.

(2) Install Anaconda according to [tutorial](http://docs.continuum.io/anaconda/install#linux-install)
```
bash ~/Downloads/Anaconda3-2.4.0-Linux-x86_64.sh
export PATH=$HOME/anaconda2/bin:$PATH
```

For the `export` line, you may want to add them to `.bashrc` or `.bash_profile` in your home directory. Otherwise you will need to type it everytime you log in.


(3) Confirm conda is installed and update conda
```
conda -V
conda update conda
```

(4) Create a virtual environment ([tutorial](http://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/)). I will call it `anaCogent`. Type `y` to agree to the interactive questions.
```
conda create -n anaCogent python=2.7 anaconda
source activate anaCogent
```

Once you have activated the virtualenv, you should see your prompt changing to something like this:

```
(anaCogent)-bash-4.1$
```

(5) Install additional libraries that we need
```
conda install -n anaCogent biopython
conda install -n anaCogent scikit-learn
```

(6) Clone ANGEL repo and install
```
git clone https://github.com/Magdoll/ANGEL.git
cd ANGEL
python setup.py build
python setup.py install
```



### Option 2: Install ANGEL using virtualenv

The alternative to Anaconda is to set up and activate a virtual environment before installation. See [here](https://github.com/PacificBiosciences/cDNA_primer/wiki/Setting-up-virtualenv-and-installing-pbtranscript-tofu) for installation details.

Once you have the virtualenv activated, install the dependencies using `pip`:

```
pip install numpy
pip install biopython
pip install scikit-learn
```


You can download this GitHub repository in many ways. Here we assume you will be using ``git clone``:

```
git clone https://github.com/Magdoll/ANGEL.git
cd ANGEL
python setup.py build
python setup.py install
```



## USAGE

<a href="dumb"/>
#### Dumb ORF prediction

`dumb_predict.py` takes as input a FASTA file. It outputs the longest ORF (or first ORF) that exceed the user-defined minimum length and have a positive log-odds scores based on hexamer frequencies.


Usage:

```
dumb_predict.py <fasta_filename> <output_prefix> 
       [--min_aa_length MIN_AA_LENGTH]
       [--use_firstORF]
       [--use_rev_strand] [--cpus CPUS]
```

See the provided training example below. By default, only the forward strand is used. This is especially true for PacBio transcriptome sequencing output that often already has correct strand.

```
cd ANGEL/training_example
dumb_predict.py test.human_1000seqs.fa test.human.dumb --min_aa_length 300 --cpus 24
```

The output consists of ``test.human.dumb.final.pep``, ``test.human.dumb.final.cds`` and ``test.human.dumb.final.utr``, which are the results of longest ORF prediction.

By default, the longest ORF (regardless of frame) is chosen as output. This longest ORF must exceed the minimum length (`--min_aa_length`) and also have a positive log-odds score. If the `--use_firstORF` option is turned on, however, then no log-odds score is computed; instead, the first ORF (regardless of frame) that exceeds the minimum length is output. The `--use_firstORF` option is useful in cases where users believe the first ORF, and not the longest ORF, is the more correct prediction.


<a name="trainset"/>
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

The output files are: ``test.human.dumb.final.training.cds``, ``test.human.dumb.final.training.utr`` and ``test.human.dumb.final.training.pep``.


<a name="training"/>
#### ANGEL classifier training


`angel_train.py` takes a CDS FASTA file and a UTR FASTA file and outputs a trained classifier pickle file. 

Usage:

```
angel_train.py <cds_filename> <utr_filename> <output_pickle> [--cpus CPUS]
```

Example:

```
angel_train.py test.human.dumb.final.training.cds test.human.dumb.final.training.utr \
      test.human.dumb.final.classifier.pickle --cpus 12
```

On a typical 500-sequence training dataset, the training may take several hours. 

**NOTE** I have found that sometimes `angel_train.py` hangs if you use more than 12 cores (regardless of memory usage), possibly due to issues with Python's multiprocessing. 


<a name="angelpredict"/>
#### ANGEL ORF prediction

Usage:

```
angel_predict.py <input_fasta> <classifier_pickle> <output_prefix>
       [--min_angel_aa_length] [--min_dumb_aa_length] 
       [--max_angel_secondORF_distance]
       [--use_rev_strand] [--output_mode [best,all]]
       [--cpus]
```

For each sequence, ANGEL uses the classifier to predict the coding potential of each codon in each of the three frames, then finds the most likely stretch of window that is the open reading frame. It then does the following:

* If there is only one predicted ORF, tag it as ``confident``. In this case, it is unlikely there are sequencing errors.
* If there are multiple ORFs but only one exceeds `min_angel_aa_length` threshold, output that one ORF and tag it as ``likely``. In this case, the program is semi-positive that there are no sequencing errors or that the error occurs at the ends of the CDS, allowing successful prediction of a relatively long continuous ORF.
* If there are multiple ORFs (which could be in multiple frames) exceeding the length threshold, output all of them (provided the later ORFs are close enough to the upstream ORFs by the threshold determined by `max_angel_secondORF_distance`) and tag them as ``suspicious``. In this case, the ORFs could be genuine complete ORFs (indicating a polycistronic transcript or alternative ORFs) or fragments of the same ORF that has frame shift due to uncorrected sequencing errors.

The longest ORF length from the ANGEL process is recorded as ``T``. Then, the same dumb ORF prediction (which simply looks for the longest stretch of ORF without stop codon interruption) is done on the forward three frames. The longest predicted dumb ORF is also outputted if, and only if, its length is greater than both ``T`` and `min_dumb_aa_length`. This is a fallback process in case the ANGEL classifier failed to detect coding potential in the CDS region.

If `--use_rev_strand` is used, then the same process is repeated on the reverse-complement of the sequence. 

The `--output_mode` option is used to determine what to output when there are ORF predictions from both ANGEL and dumb and possibly the reverse strand (if `--use_rev_strand` is on).

By default, the `--output_mode=best` will pick the longest ORF, whether that be from ANGEL or dumb, plus or reverse (if `--use_rev_strand` is on) strand.

If `--output_mode=all`, then all predictions will be output as long as they are longer than the threshold length (for ANGEL, `min_angel_aa_length` and for dumb, `min_dumb_aa_length`).

The `--max_angel_secondORF_distance` parameter is used only in the ``suspicious`` case (cases where there is more than one ORF predicted) to determine whether or not to output the downstream ORFs based on how close they are to the upstream ORFs. The reasoning is, if the ORF prediction is broken up because of sequencing errors or remaining base errors, then the two ORFs should be very close to each other (because ANGEL expects there to be only a handful of base errors, not long stretches of errors). If the second ORF is very far from the first ORF (say more than 20 bp downstream), then it is not likely to be from sequencing errors. Rather the first ORF (with its stop codon) is the correct prediction. By default `max_angel_secondORF_distance` is set to 10 bp. So any downstream ORF that is further than 10 bp away from the upstream ORF will be discarded.

Example:

```
angel_predict.py test.human_1000seqs.fa human.MCF7.random_500_for_training.pickle test.human \
      --use_rev_strand --output_mode=best --min_angel_aa_length 100 --min_dumb_aa_length 100
```

The output files are: ``<output_prefix>.ANGEL.cds``, ``<output_prefix>.ANGEL.pep``, and ``<output_prefix>.ANGEL.utr``.

Each output sequence ID has the format:
```
<seq_id> type:<tag>-<completeness> len:<ORF length (aa)> strand:<strand> pos:<CDS range>
```

Where ``tag`` is ``confident``, ``likely``, or ``suspicious`` for ANGEL predictions, and ``dumb`` for dumb ORF predictions.

``completeness`` is either ``complete``, ``5partial``, ``3partial``, or ``internal`` based on the presence or absence of start and stop codons.


<a name="errgenome"/>
#### Using genomic to error correct first before ORF prediction

If a high-quality genome is available and you want to additionally run ANGEL on a genomic version of the transcripts, you can use the `err_correct_w_genome.py` script from [cDNA_Cupcake](https://github.com/Magdoll/cDNA_Cupcake/). cDNA_Cupcake is a light-weight repository for simple manipulation scripts. See the [err_correct_w_genome tutorial](https://github.com/Magdoll/cDNA_Cupcake/wiki/Sequence-Manipulation-Wiki#errgenome) on how to obtain a genomic version first.

Once you have the genome-corrected version, you can run ANGEL separately on it using either just `dumb_predict.py` (if you trust the genome bases and alignment exon boundaries to be precise) or `angel_predict.py` (to still allow for some errors).

With two versions of ANGEL output (one from the PacBio version, one from the genomic version), you can use the script `angel_compare_files.py` to pick the longest ORF for each sequence from the two files.

```
angel_comapre_files.py <ANGEL_prefix1> <ANGEL_prefix2> <output_prefix>
```

Below is an example:

```
# run ANGEL on PacBio sequence, now obtain test1.cds, test1.utr, test1.pep
# run ANGEL on genomic sequence, now obtain test2.cds, test2.utr, test2.pep

angel_compare_files.py test1 test2 test3
```

The result will be output to: test3.cds, test3.pep, test3.utr, and test3.compare.txt


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
