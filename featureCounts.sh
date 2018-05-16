#!/bin/bash

#BSUB -J staralign
#BSUB -U rsmasw

module load subread

featureCounts -a pdam_1415_maker.gtf -o feature_counts.out -t exon -g gene_id bamlist.txt
