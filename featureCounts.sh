#!/bin/bash

#PBS -d Set_working_directory_to_where_bam_files_are_stored
#PBS -lnodes=1:ppn=4,pmem=8gb
#PBS -lwalltime=48:00:00

genfilework=Homo_sapiens.GRCh38.dna.primary_assembly.reordered.fa
gtffile=Homo_sapiens.GRCh38.gtf
bam_file_to_count=bam_files_to_count.txt # contains complete paths to bam files in one line, ex:
# bamfile1.bam bamfile2.bam bamfile3.bam

featureCounts -T 4 -g gene_id -J -G $genfilework -p -s 2 -d 20 -a $gtffile -o counts.txt $(eval 'cat $bam_file_to_count')
