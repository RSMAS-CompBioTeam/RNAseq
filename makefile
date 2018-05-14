#run DESeq on featurecounts output
Rscript Basic_SNP_analyses.R

#install samtools
wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
tar -xvjf samtools-1.8.tar.bz2
cd samtools-1.8
make
cd ..

#download fasta of protein models
wget TBD

#use samtools to extract DE genes
xargs samtools-1.8/samtools faidx pdam_1415_maker.faa < DE_genes.txt > DE_genes.faa
