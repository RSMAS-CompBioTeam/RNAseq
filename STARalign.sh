#!/bin/bash
#~/scripts/EAPSI.WT2-master/STARalign.sh
#/scratch/projects/transcriptomics/mikeconnelly/scripts/EAPSI.WT-master/STARalign.sh
#purpose: align trimmed RNAseq reads using STAR on the Pegasus bigmem queue

#BSUB -J staralign
#BSUB -q bigmem
#BSUB -P transcriptomics
#BSUB -o star%J.out
#BSUB -e star%J.err
#BSUB -n 8
#BSUB -B m.connelly@rsmas.miami.edu
#BSUB -N

#specify variable containing sequence file prefixes and directory paths
mcs="/scratch/projects/transcriptomics/mikeconnelly"
coldir="/scratch/projects/transcriptomics/mikeconnelly/sequences/EAPSI/wanglitung"
exp="alltreatments"
EAPSIsamples="Wt1-1a Wt1-1b Wt1-1c Wt1-2a Wt1-2b Wt1-3a Wt1-3b Wt1-3c Wt1-4a Wt1-4b Wt1-4c Wt1-5a Wt1-5b Wt1-5c Wt1-6b Wt1-6c Wt2-1a Wt2-1b Wt2-1c Wt2-2b Wt2-2c Wt2-3a Wt2-3b Wt2-3c Wt2-4a Wt2-5a Wt2-5b Wt2-5c Wt2-6a Wt2-6b Wt2-6c"

#Run STAR aligner
echo "These are the reads to be aligned to the Pocillopora reference genome: $EAPSIsamples"
for EAPSIsample in $EAPSIsamples
do \
echo "Aligning ${EAPSIsample}" 
${mcs}/programs/STAR-2.5.3a/bin/Linux_x86_64/STAR \
--runMode alignReads \
--runThreadN 16 \
--readFilesIn ${coldir}/${exp}/trimmomaticreads/${EAPSIsample}_trim.fastq.gz \
--readFilesCommand gunzip -c \
--genomeDir ${mcs}/sequences/genomes/coral/pocillopora/STARindex \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile  ${mcs}/sequences/genomes/coral/pocillopora/pdam_genome.gff \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${coldir}/${exp}/STARalign/${EAPSIsample} 
#lets me know file is done
echo "STAR alignment of $EAPSIsample complete"
done

#Call next scripts in analysis pipeline
bsub -P transcriptomics < /scratch/projects/transcriptomics/mikeconnelly/scripts/EAPSI.WT-master/STARfC.sh