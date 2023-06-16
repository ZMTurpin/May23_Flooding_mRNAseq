#!/bin/bash

#trim contaminating adapter sequences with cutadapt
        #usage cutadapt -a <adapterFwd> -A<adapterRev> -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

        #create I/O filename arrays
#mkdir trim

R1=($(ls Maize_Flooding_RNA-seq_pool/*R1*))
R2=($(ls Maize_Flooding_RNA-seq_pool/*R2*))

R1out=($(echo ${R1[@]} | awk '{gsub("R1_001","R1_clip") ; print $0 }'))

R2out=($(echo ${R2[@]} | awk '{gsub("R2_001","R2_clip") ; print $0 }'))

#AdaptorRead 1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#AdaptorRead 2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

for (( a=0 ; a<${#R1[@]} ; a++))
        do
        cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o ${R1out[$a]} -p ${R2out[$a]} ${R1[$a]} ${R2[$a]}
        done

#run fastqc on trimmed reads

trim=($(for a in $(ls trim/*clip*) ; do printf $a" " ; done))

mkdir ClipFASTQC
fastqc -t 6 -o ClipFASTQC ${trim[@]}



