#run DESeq on featurecounts output
Rscript Basic_SNP_analyses.R

#download fasta of protein models
wget TBD

#use samtools to extract DE genes
xargs samtools-1.8/samtools faidx pdam_1415_maker.faa < DE_genes.txt > DE_genes.faa
