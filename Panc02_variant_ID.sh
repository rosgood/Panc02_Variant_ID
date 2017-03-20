#!/bin/bash

#to make bowtie indeces:
bowtie2-build mm9.fa mm9_bowtie2

mkdir tmp/

#bowtie2 --fast -x mm9_bowtie2 -q -1 s_6_1.fastq -2 s_6_2.fastq | samtools_0.1.18 view - -Sb -h -T mm9.fa -o tmp/bt2_s_6.bam;

#pipe directly to samtools
bowtie2 -x mm9_bowtie2 -q -1 sorted_s_6_1.fq -2 sorted_s_6_2.fq | samtools_0.1.18 view -bS - > tmp/bt2_s_6.bam

#after align by bowtie2, convert to sorted, indexed bams, then do calls

cd tmp

#samtools_0.1.18 sort bt2_s_6.bam sort_bt2_s_6;
#samtools_0.1.18 index sort_bt2_s_6.bam;

#coord sort and add RGs for potential MuTect calling

java -jar /home/heatherkinkead/picard/SortSam.jar I=bt2_s_6.bam O=coord_sort_bt2_s_6.bam SO=coordinate;
java -jar /home/heatherkinkead/picard/AddOrReplaceReadGroups.jar I=coord_sort_bt2_s_6.bam O=RG_sort_bt2_s_6.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1;

samtools_0.1.18 index RG_sort_bt2_s_6.bam;
#Having problems with freebays reading bam file, so convert to sam, cat header to sam, then convert back to back. samtools reheader not working
#samtools_0.1.18 view coord_sort_bt2_s_6.bam > test.sam
#cat header tmp/test.sam > test_w_header.sam
#samtools_0.1.18 view -bSh test_w_header.sam > tmp/bt2_s_6_w_head.bam

freebayes -f ~/data/Panc02_exome/mm9.fa --genotype-qualities ~/data/Panc02_exome/tmp/RG_sort_bt2_s_6.bam | vcfqualfilter -c 20 > bt2_Panc02_filtered.vcf;

#now follow annotation pipeline

files=*_filtered.vcf

for i in $files
do
	convert2annovar.pl --format vcf4old --includeinfo --comment "$i" > "${i%_filtered.vcf}"_av_ann;
	annotate_variation.pl --geneanno --buildver mm9 --comment --separate --webfrom ucsc "${i%_filtered.vcf}"_av_ann ~/annovar/mousedb/
done;

#get only nonsyn. vars
files=*.exonic_variant_function
for i in $files
do
	grep "nonsynonymous SNV" "$i" > "${i%.exonic_variant_function}"_missense.exonic_variant_function
#to get annotation statistics
	cut -f 2 "$i" | sort | uniq -c | sort -n > "${i%_missense.exonic_variant_function}"_exonic_variant_stats.txt
done;

#convert exonic variants file to amino acid fasta
files=*_missense.exonic_variant_function
for i in $files
do
###potential fasta errors leading to downstream errors, be sure to check all muts in IGV at end
	coding_change.pl --includesnp "$i" /home/heatherkinkead/annovar/mousedb/mm9_refGene.txt /home/heatherkinkead/annovar/mousedb/mm9_refGeneMrna.fa > "${i%_missense.exonic_variant_function}"_coding_seqs.fasta
#get mutation matched to unique ID
        #this potentially leading to multiple lines from multiple overlapping genes
	#cut -f 3 "$i" | tr "," "\n" | grep -v "^$" | tr ":" "\t" > "${i%_missense.exonic_variant_function}"_matched_names_ref_gene.txt
	cut -f 3 "$i" | sed 's/,.*//g' | grep -v "^$" | uniq | tr ":" "\t" > "${i%_missense.exonic_variant_function}"_matched_names_ref_gene.txt
done;

#create a tab delimited file to get rid of wt seqs
files=*_coding_seqs.fasta
for i in $files
do
##There is a problem with this awk command
        awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' "$i" > "${i%_av_ann_coding_seqs.fasta}"_table
#create headers for ID creation
        grep ">" "$i" | grep -v "WILDTYPE" > "${i%.fasta}"_mut_headers.txt
done;

#get rid of wt seqs
files=*_table
for i in $files
do
        grep -v "WILDTYPE" "$i" > "${i%_table}"_muts_prots_only
done;

#get rid of first column
files=*_muts_prots_only
for i in $files
do
#in case of duplicate calls which screw up the R table, add in a uniq step
        #cut -d ' ' -f 2-12 "$i" > "$i"_fixed
	cut -d ' ' -f 2-12 "$i" | uniq >"$i"_fixed
done;

#create unique ID
files=*_headers.txt
for i in $files
do
    cut -d ' ' -f 2-3 "$i" > "${i%_mutect_overlap_av_ann_coding_seqs_mut_headers.txt}"_unique_ID.txt
done;


#don't think I need to sort both lists so that they can be joined in R
#open R and begin merging tables
R CMD BATCH "/home/heatherkinkead/data/Panc02_exome/R_merge_tables.R";
# test_R_script.R merges tables together, once that is done pick up here

#get mut locs
files=*_best_ID_R
for i in $files
do
        cut -f 3 "$i" | sed 's/[^0-9]//g' | tail -n +1 > "${i%_best_ID_R}"_mut_locs
	cut -f 1 "$i" | tail -n +1 > "${i%_best_ID_R}"_mut_locs_ID
        paste "${i%_best_ID_R}"_mut_locs_ID "${i%_best_ID_R}"_mut_locs > "${i%_best_ID_R}"_mut_locs_w_ID
done;

#create fasta file
files=*_best_ID_R
for i in $files
do
        awk 'BEGIN {FS="\t"} {print ">"$1"\n"$4}' "$i" | fold -w 80 > "$i".fasta
done;

#get peptides with R
#need to do this individually until can figure out code
#Use R_regions_peptides_indivly.R
R CMD BATCH "/home/heatherkinkead/data/Panc02_exome/R_regions_peptides_indivly.R"

grep -A 2 Warning "/home/heatherkinkead/data/Panc02_exome/R_regions_peptides_indivly.Rout"

#cat files together from separate reads only if present in both
#uniq -d to get only in both files
##do netMHC calls with Exome_NetMHC_predictions

touch ~/data/Panc02_exome/tmp/bt2_Panc02_peptides_NetMHC.out

files=*_peptides.fasta
for i in $files
do
	netMHC -a H-2-Db,H-2-Kb -l 8 -s ~/data/Panc02_exome/tmp/bt2_Panc02_peptides.fasta >> bt2_Panc02_peptides_NetMHC.out
	netMHC -a H-2-Db,H-2-Kb -l 9 -s ~/data/Panc02_exome/tmp/bt2_Panc02_peptides.fasta >> bt2_Panc02_peptides_NetMHC.out
	netMHC -a H-2-Db,H-2-Kb -l 10 -s ~/data/Panc02_exome/tmp/bt2_Panc02_peptides.fasta >> bt2_Panc02_peptides_NetMHC.out
done

#get rid of awful formatting
#grep 'NM_' T1A_vs_T4_results.out | sort -k4n > T1AT4.clean
grep H-2 bt2_Panc02_peptides_NetMHC.out | grep -v version | sort -k4n > bt2_Panc02_peptides.clean

#get rid of WB awful column
#sed 's/WB//' T1AT4.clean | sed 's/SB//' > T1AT4.reallyclean
sed 's/WB//' bt2_Panc02_peptides.clean | sed 's/SB//' > bt2_Panc02_peptides.reallyclean
#Use R to clean files
R CMD BATCH "/home/heatherkinkead/data/Panc02_exome/NetMHC_clean.R"

#this command doesn`t have affinity sorted properly, but won`t miss any epitopes
#this was how the epitopes were run to generate the file, so they don't have the best predicted value. filled in by hand later
#sort -u -k5,5 bt2_Panc02_peptides_best.predicted > bt2_uniq_Panc02_predictions
sort -k5,5 bt2_Panc02_peptides_best.predicted | sort -u -k4,4n > bt2_uniq_Panc02_predictions
#cat header on so R can merge columns
cat ~/data/Panc02_exome/NetMHC_header bt2_uniq_Panc02_predictions > bt2_uniq_predicts_w_header

#add peptides onto predictions
R CMD BATCH "/home/heatherkinkead/data/Panc02_exome/R_get_full_peps_predicts.R"

