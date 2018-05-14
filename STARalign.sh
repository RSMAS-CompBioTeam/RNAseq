#!/bin/bash
#adapted from a script by Mike Connelly
#purpose: align trimmed RNAseq reads using STAR on the Pegasus bigmem queue

#BSUB -J staralign
#BSUB -U rsmasw
#BSUB -o star%J.out
#BSUB -e star%J.err

module load star

for i in $@
do \
BASE=$(basename i .fq.gz)
echo "Aligning $i" 
STAR \
--runMode alignReads \
--readFilesIn $i \
--readFilesCommand gunzip -c \
--genomeDir STARindex \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile pdam_1415_maker.gtf \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $BASE
#lets me know file is done
echo "STAR alignment of $i complete"
done
