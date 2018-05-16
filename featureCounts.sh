#!/bin/bash

#BSUB -J staralign
#BSUB -U rsmasw
#BSUB -o star%J.out
#BSUB -e star%J.err

module load subread

featureCounts -a pdam_1415_maker.gtf -o feature_counts.out -t exon -g gene_id $@
