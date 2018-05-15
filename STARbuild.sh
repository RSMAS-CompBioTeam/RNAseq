#!/bin/bash
#~/RSMAS-CompBioTeam/RNAseq/STARbuild.sh
#purpose: build genome indices for use with the STAR aligner

#BSUB -J starbuild
#BSUB -U rsmasw

module load star

mkdir STARindex

STAR \
--runMode genomeBuild \
--genomeDir STARindex \
--genomeFastaFiles pdam_genome_1415.fasta \
--sjdbGTFfile pdam_1415_maker.gtf \
