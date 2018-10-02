#!/bin/bash

#PBS -d Path_containing_fastq_files_with_SRA_accession_PRJNA491629
#PBS -lnodes=1:ppn=1,pmem=32gb
#PBS -lwalltime=2:00:00:00
#PBS -t 1-5
i=$PBS_ARRAYID

genfilework=Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtffile=Homo_sapiens.GRCh38.gtf
genomeloc=path_to_STAR_index

lane=L2_INDEX
read1=140908_I312_FCC56NNACXX_$lane$i\_1.fq
read2=140908_I312_FCC56NNACXX_$lane$i\_2.fq

STAR --clip3pNbases 2 --clip5pNbases 13 --outFileNamePrefix ./$lane$i --genomeDir $genomeloc --sjdbGTFfile $gtffile --outSAMtype BAM Unsorted SortedByCoordinate --readFilesCommand zcat --readFilesIn $read1 $read2 --runThreadN 1
samtools index $lane$i\Aligned.sortedByCoord.out.bam

# same for files with L7 index
lane=L7_INDEX
read1=140908_I312_FCC56NNACXX_$lane$i\_1.fq
read2=140908_I312_FCC56NNACXX_$lane$i\_2.fq

STAR --clip3pNbases 2 --clip5pNbases 13 --outFileNamePrefix ./$lane$i --genomeDir $genomeloc --sjdbGTFfile $gtffile --outSAMtype BAM Unsorted SortedByCoordinate --readFilesCommand zcat --readFilesIn $read1 $read2 --runThreadN 1
samtools index $lane$i\Aligned.sortedByCoord.out.bam
