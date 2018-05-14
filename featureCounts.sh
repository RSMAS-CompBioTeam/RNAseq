#!/bin/bash

featureCounts -a pdam_1415_maker.gtf -o feature_counts.out -t exon -g gene_id $@
