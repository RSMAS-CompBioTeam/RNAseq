#!/bin/bash
#adapted from a script by Mike Connelly
#purpose: align trimmed RNAseq reads using STAR on the Pegasus bigmem queue

#BSUB -J staralign
#BSUB -U rsmasw

module load star

while read i; do
BASE=$(basename $i .fq.gz)
echo "Aligning $i" 
STAR \
--runMode alignReads \
--readFilesIn $i \
--readFilesCommand gunzip -c \
--genomeDir STARindex \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile pdam_1415_maker.gtf \
--outFileNamePrefix "$BASE"
#lets me know file is done
echo "STAR alignment of $i complete"
#index the bam file
samtools view -b "$BASE"Aligned.out.sam > "$BASE"_unsorted.bam
samtools sort "$BASE"_unsorted.bam > "$BASE".bam
samtools index "$BASE".bam
rm "$BASE"Aligned.out.sam
rm "$BASE"_unsorted.bam
done < samples.txt
