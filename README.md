ANGEL: Robust Open Reading Frame prediction
=====

The program is divided into three parts:

-- dumb ORF prediction: ORF prediction that outputs all longest ORFs in all frames. Can be used to create a top training dataset.

-- ANGEL classifier training: train a coding potential classifier based on given training data. Must provide both CDS and UTR for positive and negative training set. 

-- ANGEL robust ORF prediction: ORF prediction based on both the ANGEL classifier and the dumb ORF prediction.


The ANGEL classifer is based on the classifer described in the Shimizu *et al.*, "**ANGLE: a sequencing errors resistant program for predicting protein coding regions in unfinished cDNA.**", *J Bioinform Comput Biol*, (2006). The naming change to ANGEL is intentional to disambiguate from the author implementation of ANGLE.



