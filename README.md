# SA GWAS

#### 1. Phenotype and quant. gen analysis

-Male and female fitness data QC

-Male and female heritability and rmf estimation using MCMCglmm

-SNP heritability of the antagonism index

-Define rotated matrix to map antagonism/concordant fitness effects of each line

-Male and female fitness plot + antagonism index overlay

#### 2. Genotype quality control

-Depth and GQ filtering

-MAF and call rate filtering

-PCA including outgroup population (DGRP) to identify outliers

#### 3. LHm pop gen

-Examination of pairwise relatedness between individuals

-LD in LHm

#### 4. Association mapping

-LDAK mixed model association analysis

-Manhattan plot

-QQ-plot

-P-value histogram

-FDR estimation

-Comparison of antagonistic and concordant P-values and Q-values

-Permutation-based analysis

-LD clumping to identify number of independent hits

#### 5. Functional analyses

-Partitioning heritability into autosome vs. X

-Antagonistic vs. non-antagonistic candidates on autosome vs. X 

-Partitioning heritability into variant effect categories

-Antagonistic vs. non-antagonistic candidates among variant effect categories

-Antagonistic genes: GO analysis

-Antagonistic genes: sex-biased expression

-Antagonistic genes: pleiotropy

#### 6. SNP-based bal sel analyses (Dmel)

-Convert Genome Nexus dasta files to vcf format and r6 coordinates

-Estimate minor allele frequencies in three populations (DGRP, ZI, SA), including sites that are monomorphic in the population of interest but polymorphic in LHm (MAF=0)

-Incoporate linked selection / recombination rate estimates

-Plot MAF spectrum at LD-pruned antagonistic vs. non-antagonistic SNPs

-Comparison of antagonistic and non-antagonistic MAF while correcting for MAF ascertainemnt bias and linked selection ('Analysis A')

-Relationship between GWAS effect size and probability that SNP is polymorphic, correcting for MAF ascertainemnt bias and linked selection ('Analysis B')

-Relationship between GWAS effect size and MAF, correcting for MAF ascertainemnt bias and linked selection ('Analysis C')

#### 7. Window- and LD-based bal sel analyses

-Import 'antagonistic windows'

-Calculate Tajima's D values in sliding windows

-Calculate Fst values in sliding windows

-Model Tajima's D and Fst as a function of antagonistic/non-antagonistic status + include LS as covariate

-Estimate LD between pairs of sites in ZI 

-Compare pairwise LD between antagonistic sites relative to non-antagonistic sites

#### 8. Data files

-Code used to produce summary data files. These data files are sufficient to reproduce the figures presented in the manuscript.

