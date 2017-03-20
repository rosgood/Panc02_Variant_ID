# Annotating Variants from Whole Exome Sequencing Data

This pipeline assumes bedtools, samtools, picardtools, freebayes, annovar, and NetMHC3.4 are all installed as local packages
Reference genome file locations need to be updated in the script to reflect correct locations

## Starting Files

Sorted .bam files first need to be converted to .fastq for bowtie2 alignment

`sh bam2fastq.sh`

## Variant Annotation and Peptide Sequence Extraction and Prediction

Verify all pathnames

*`Panc02_variant_ID.sh`
*`R_merge_tables.R`
*`R_get_full_peps_predicts.R`
*`R_regions_peptides_indivly.R`
*`NetMHC_clean.R`


Then run
`sh Panc02_variant_ID.sh`
