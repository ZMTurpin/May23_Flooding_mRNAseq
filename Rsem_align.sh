#!/bin/bash

#Create STAR indices and align to transcript models in RSEM. Estimate transcript- and gene- level mRNA abundances


#align Flooding RNA-seq reads to transcriptome and quantify TPM abundance for each known transcript model


#prepare STAR references with RSEM (set sjdboverhang = 49)

#Example usage for 1 assembly
#rsem-prepare-reference -p 5 --star --star-sjdboverhang 49 --gff3 /psi/basslab/zturpin/genomes/Zea_mays/B73_v5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 /psi/basslab/zturpin/genomes/Zea_mays/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa ./Zm-B73v5


gffs=($(ls /psi/basslab/zturpin/genomes/Zea_mays/*/*gff3))
refs=($(ls /psi/basslab/zturpin/genomes/Zea_mays/*/*fa | grep -v "chr"))
names=($(printf %s"\n" ${gffs[@]} | awk '{split($0,a,"/") ; print a[7]}'))

for ((a=0 ; a<${#gffs[@]} ; a++ ))
do
	mkdir ${names[$a]}_rsemReference
	rsem-prepare-reference -p 5 --star --star-sjdboverhang 49 --gff3 ${gffs[$a]} ${refs[$a]} ${names[$a]}_rsemReference/${names[$a]}
done

#calculate expression (RSEM)
#NOTE RSEM HAS A BUG! IT DOES NOT SUPPORT ON-THE-FLY UNZIPPING AND PASSING ON STDOUT TO STAR. workaround by unzipping all fastqs and deleting unzipped files after alignment
#Example usage for 1 replicate

#rsem-calculate-expression --strandedness reverse -p 5 --star --append-names --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files --paired-end /psi/bassdata/zturpin/Flooding/trim/B73_control_1_S24_L002_R1_clip.fastq /psi/bassdata/zturpin/Flooding/trim/B73_control_1_S24_L002_R2_clip.fastq ./Zm-B73v5 B73_control_1

gunzip -k /psi/bassdata/zturpin/Flooding_MOAseq/cat_trim/*gz

R1=($(ls /psi/bassdata/zturpin/Flooding_MOAseq/cat_trim/*R1_trim.fastq.gz))
R2=($(ls /psi/bassdata/zturpin/Flooding_MOAseq/cat_trim/*R2_trim.fastq.gz))

sampleNames=($(printf %s"\n" ${R1[@]} | awk '{split($0,a,"/") ; split(a[7],b,"_") ; print b[1]"_"b[2]"_"b[3]}'))
#refArray=($(for b in ${names[@]} ; do printf %s"\n" $b_refArray $b $b $b $b $b ; done ))
refArray=($(for b in ${names[@]} ; do printf %s"\n" $b"_rsemReference/"$b $b"_rsemReference/"$b $b"_rsemReference/"$b $b"_rsemReference/"$b $b"_rsemReference/"$b $b"_rsemReference/"$b ; done))
for ((c=0; c<${#R1[@]}; c++))
do
	rsem-calculate-expression --strandedness reverse -p 5 --star --append-names --output-genome-bam --sort-bam-by-coordinate --keep-intermediate-files --paired-end ${R1[$c]} ${R2[$c]} ${refArray[$c]} ${sampleNames[$c]}

done

unalias rm
rm /psi/bassdata/zturpin/Flooding/trim/*fastq
