# pADamID_ABCB1_Novaseq
This pipeline has been used to process pA-DamID samples sequenced with Novaseq platform.

pA-DamID scripts
All the scripts for pA-DamID data processing and analyses for Novaseq platform.

The Snakemake pipeline can be found in "bin", where "snakemake/damid.snake" contains the actual snakemake pipeline compatibe with fastq files geenrated with Novaseq and "snakemake/config_ste_20220704Anna.yaml" contains the sample information and parameters. Further processing of the data was done with R markdown files:

SGM20220725_Anna_Analysis_TAxol_Resistance_Novaseq.Rmd: this script loads data in R and assign a Nl association score to genes. This script was used to generate  scatter plots for genes in many different conditions. Finally, the script was used to calculate correlations between independent replicates and significant detachment on ABCB1 gene.
