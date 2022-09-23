# Genome-wide sexually antagonistic variants reveal long-standing constraints on sexual dimorphism in fruit flies


Welcome!

This github page contains the code underlying analyses and figures presented in 'Genome-wide sexually antagonistic variants reveal long-standing constraints on sexual dimorphism in fruit flies' (Ruzicka et al. 2019, PLOS Biol). 

The code is organised into text files 1-9, which approximately follow the order of the analyses presented in the manuscript. The text files are usually a mixture of R and UNIX scripts. The contents of each file are described in bullet-point form below. We have also annotated the code for clarity, but if it's insufficiently clear please e-mail me at filip.ruzicka [at] monash.edu.

The 9th text file contains more detailed descriptions of the datafiles shared on zenodo (https://zenodo.org/record/2623225). These zenodo datafiles contain the raw data used to produce the figures in the manuscript, and will be the most useful objects for future users.


#### 1. Phenotype and quant. gen analysis

-Male and female fitness data QC

-Male and female heritability and rmf estimation using MCMCglmm

-Define rotated matrix to map antagonism/concordant fitness effects of each line

-Male and female fitness plot + antagonism index overlay

#### 2. Genotype quality control

-Depth and GQ filtering

-MAF and call rate filtering

-PCA including outgroup population (DGRP) to identify outliers

#### 3. LHm pop gen

-Pairwise relatedness between individuals

-LD in LHm

#### 4. Association mapping

-SNP heritability of the antagonism index

-LDAK mixed model association analysis

-Manhattan plot

-QQ-plot

-FDR estimation

-Permutation-based analysis

-LD clumping to identify number of independent hits

-Comparison of antagonistic and concordant P-values and Q-values

-Degree of clustering of antagonistic SNPs

#### 5. Functional analyses

-Partitioning heritability into autosome vs. X

-Partitioning heritability into variant effect categories

-Antagonistic vs. non-antagonistic candidates on autosome vs. X 

-Antagonistic vs. non-antagonistic candidates among variant effect categories

-Antagonistic genes: sex-biased expression

-Antagonistic genes: pleiotropy

#### 6. SNP-based bal sel analyses (Dmel)

-Convert Genome Nexus fasta files to vcf format and r6 coordinates

-Estimate MAFs in three populations (DGRP, ZI, SA)

-Incoporate linked selection / recombination rate estimates

-Plot MAF spectrum at LD-pruned antagonistic vs. non-antagonistic SNPs

-Comparison of antagonistic and non-antagonistic MAF while correcting for MAF ascertainemnt bias and linked selection ('Analysis A')

-Relationship between GWAS effect size and probability that SNP is polymorphic, correcting for MAF ascertainemnt bias and linked selection ('Analysis B')

-Relationship between GWAS effect size and MAF, correcting for MAF ascertainemnt bias and linked selection ('Analysis C')

#### 7. Window- and LD-based bal sel analyses

-Conduct window-based association test

-Calculate Tajima's D values in sliding windows

-Calculate Fst values in sliding windows

-Model Tajima's D and Fst as a function of antagonistic/non-antagonistic status + include LS as covariate

-Estimate LD between pairs of sites in ZI 

-Compare pairwise LD between antagonistic and non-antagonistic sites in ZI

#### 8. SNP-based bal sel analyses (Dsim, Dyak)

-Import allele frequencies for two D. simulans datasets and one D. yakuba dataset

-Relationship between GWAS effect size and probability that SNP is polymorphic, correcting for MAF ascertainemnt bias and linked selection ('Analysis B')

-Comparison of antagonistic and non-antagonistic trans-specific status while correcting for MAF ascertainemnt bias and linked selection ('Analysis A')

#### 9. Data files

-Description of summary data files 

-Code used to produce summary data files. These data files should be sufficient to reproduce the figures presented in the manuscript

