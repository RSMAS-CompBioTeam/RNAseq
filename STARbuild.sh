#!/bin/bash
#~/RSMAS-CompBioTeam/RNAseq/STARbuild.sh
#purpose: build genome indices for use with the STAR aligner

#BSUB -J starbuild
#BSUB -U rsmasw
#BSUB -P ccsfellows

module load star

mkdir STARindex

STAR \
--runMode genomeGenerate \
--genomeDir STARindex \
--genomeFastaFiles pdam_genome_1415.fasta \
--sjdbGTFfile pdam_1415_maker.gtf
