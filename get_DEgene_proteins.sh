#!/bin/bash
#USAGE: bash get_DEgene_proteins.sh DEgenes.txt proteins.faa output.faa

xargs samtools-1.8/samtools faidx $2 < $1 > $3
