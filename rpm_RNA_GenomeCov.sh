#!/bin/bash
#This script prepares rpm-normalized per-base coverage files (bw) for each replicate of an RNA-seq experiment AND merges the replicates.
#inputs are STAR-aligned read pairs
#ZMT_JULY27_2023

#Create arrays of genome.sorted.bam file paths (from RSEM)

B73bams=($(ls B73*genome.sorted.bam))
Ohbams=($(ls Oh*genome.sorted.bam))
Kybams=($(ls Ky*genome.sorted.bam))
Sbbams=($(ls SB*genome.sorted.bam))

B73root=($(echo ${B73bams[@]} | awk '{gsub(".genome.sorted.bam","") ; print $0}'))
Ky21root=($(echo ${Kybams[@]} | awk '{gsub(".genome.sorted.bam","") ; print $0}'))
OhRoot=($(echo ${Ohbams[@]} | awk '{gsub(".genome.sorted.bam","") ; print $0}'))
Sbroot=($(echo ${Sbbams[@]} | awk '{gsub(".genome.sorted.bam","") ; print $0}'))

#Paths to chromsize files
B73=/psi/basslab/zturpin/genomes/Zea_mays/B73_v5/B73v5_chromsizes_chr
Ky21=/psi/basslab/zturpin/genomes/Zea_mays/Ky21/Zm-Ky21-REFERENCE-NAM-1.0_chromsizes_chr.txt
Oh7b=/psi/basslab/zturpin/genomes/Zea_mays/Oh7b/Zm-Oh7B-REFERENCE-NAM-1.0_chromsizes_chr.txt
Sorg=/psi/basslab/zturpin/genomes/Zea_mays/SorBi_NCBIv3.55/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_chromsizes_chr.txt


#B73
#compute normalization factor for each replicate
B73factors=($(for i in ${B73bams[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))

#compute normalized genomeCoverage for replicates
for (( a=0 ; a<${#B73factors[@]} ; a++))
do
	bedtools genomecov -bg -scale ${B73factors[$a]} -split -ibam ${B73bams[$a]} | grep -v "scaf" >> ${B73root[$a]}_chr_rpm.bg
done

#convert bedgraphs to bigwigs

B73bgs=($(ls B73*bg))

for (( b=0 ; b<${#B73bgs[@]} ; b++))
do
	bedGraphToBigWig ${B73bgs[$b]} $B73 ${B73root[$b]}_chr_rpm.bw
done


#Ky21
Kyfactors=($(for i in ${Kybams[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))

for (( a=0 ; a<${#Kyfactors[@]} ; a++))
do
	bedtools genomecov -bg -scale ${Kyfactors[$a]} -split -ibam ${Kybams[$a]} | grep -v "scaf" >> ${Ky21root[$a]}_chr_rpm.bg
done

Kybgs=($(ls Ky*bg))

for ((b=0 ; b<${#Kybgs[@]} ; b++))
do
	bedGraphToBigWig ${Kybgs[$b]} $Ky21 ${Ky21root[$b]}_chr_rpm.bw
done

#Oh7b

Ohfactors=($(for i in ${Ohbams[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))

for (( a=0 ; a<${#Ohfactors[@]} ; a++))
do
	bedtools genomecov -bg -scale ${Ohfactors[$a]} -split -ibam ${Ohbams[$a]} | grep -v "scaf" >> ${OhRoot[$a]}_chr_rpm.bg
done

Ohbgs=($(ls Oh*bg))
for ((b=0 ; b<${#Ohbgs[@]} ; b++))
do
	bedGraphToBigWig ${Ohbgs[$b]} $Oh7b ${OhRoot[$b]}_chr_rpm.bw
done

#Sorghum

Sbfactors=($(for i in ${Sbbams[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))

for (( a=0 ; a<${#Sbfactors[@]} ; a++))
do
        bedtools genomecov -bg -scale ${Sbfactors[$a]} -split -ibam ${Sbbams[$a]} | grep -v "super" >> ${Sbroot[$a]}"_chr_rpm.bg"
done


Sbbgs=($(ls SB*bg))

for ((b=0 ; b<${#Sbbgs[@]} ; b++))
do
        bedGraphToBigWig ${Sbbgs[$b]} $Sorg ${Sbroot[$b]}"_chr_rpm.bw"
done

#merge replicates, normalize, and compute coverage

B73c=($(ls B73_control*genome.sorted.bam))
B73f=($(ls B73_flood*genome.sorted.bam))
Ohc=($(ls Oh_control*genome.sorted.bam))
Ohf=($(ls Oh_flood*genome.sorted.bam))
Kyc=($(ls Ky_control*genome.sorted.bam))
Kyf=($(ls Ky_flood*genome.sorted.bam))
Sbc=($(ls SB_control*genome.sorted.bam))
Sbf=($(ls SB_flood*genome.sorted.bam))

samtools merge -o B73_control_merge.bam ${B73c[@]}
samtools merge -o B73_flood_merge.bam ${B73f[@]}

samtools merge -o Oh7b_control_merge.bam ${Ohc[@]}
samtools merge -o Oh7b_flood_merge.bam ${Ohf[@]}

samtools merge -o Ky21_control_merge.bam ${Kyc[@]}
samtools merge -o Ky21_flood_merge.bam ${Kyf[@]}

samtools merge -o SB_control_merge.bam ${Sbc[@]}
samtools merge -o SB_flood_merge.bam ${Sbf[@]}

#calculate merged scaling factors
B73m=($(ls B73*merge.bam))
Ohm=($(ls Oh*merge.bam))
Kym=($(ls Ky*merge.bam))
Sbm=($(ls SB*merge.bam))

Bout=($(echo ${B73m[@]} | awk '{gsub(".bam",""); print $0}'))
Oout=($(echo ${Ohm[@]} | awk '{gsub(".bam",""); print $0}'))
Kout=($(echo ${Kym[@]} | awk '{gsub(".bam",""); print $0}'))
Sout=($(echo ${Sbm[@]} | awk '{gsub(".bam",""); print $0}'))

B73Mfactors=($(for i in ${B73m[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))
OhMfactors=($(for i in ${Ohm[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))
KyMfactors=($(for i in ${Kym[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))
SbMfactors=($(for i in ${Sbm[@]} ; do samtools flagstat $i | grep "with itself and mate mapped" | awk '{print 1000000/$1}' ; done))

#compute coverage for merged and convert bedgraphs to bigwigs

#B73 merged
for ((a=0 ; a<${#B73Mfactors[@]} ; a++))
do
	bedtools genomecov -bg -scale ${B73Mfactors[$a]} -split -ibam ${B73m[$a]} | grep -v "scaf" >> ${Bout[$a]}_merge_chr_rpm.bg
done

B73mbgs=($(ls B73*merge*bg))

for ((b=0 ; b<${#B73mbgs[@]} ; b++))
do
	bedGraphToBigWig ${B73mbgs[$b]} $B73 ${Bout[$b]}_merge_chr_rpm.bw
done

#Oh7b merged
for (( a=0 ; a<${#OhMfactors[@]} ; a++))
do
        bedtools genomecov -bg -scale ${OhMfactors[$a]} -split -ibam ${Ohm[$a]} | grep -v "scaf" >> ${Oout[$a]}_merge_chr_rpm.bg
done

Ohmbgs=($(ls Oh*merge*bg))
for ((b=0 ; b<${#Ohmbgs[@]} ; b++))
do
        bedGraphToBigWig ${Ohmbgs[$b]} $Oh7b ${Oout[$b]}_merge_chr_rpm.bw
done

#Ky merged
for ((a=0 ; a<${#KyMfactors[@]} ; a++))
do
        bedtools genomecov -bg -scale ${KyMfactors[$a]} -split -ibam ${Kym[$a]} | grep -v "scaf" >> ${Kout[$a]}_merge_chr_rpm.bg
done

Ky21mbgs=($(ls Ky*merge*bg))

for ((b=0 ; b<${#Ky21mbgs[@]} ; b++))
do
        bedGraphToBigWig ${Ky21mbgs[$b]} $Ky21 ${Kout[$b]}_merge_chr_rpm.bw
done


#Sorghum Bicolor merged


for (( a=0 ; a<${#SbMfactors[@]} ; a++))
do
        bedtools genomecov -bg -scale ${SbMfactors[$a]} -split -ibam ${Sbm[$a]} | grep -v "super" >> ${Sout[$a]}"_merge_chr_rpm.bg"
done


Sbmbgs=($(ls SB*merge*bg))

for ((b=0 ; b<${#Sbmbgs[@]} ; b++))
do
        bedGraphToBigWig ${Sbmbgs[$b]} $Sorg ${Sout[$b]}"_chr_rpm.bw"
done


