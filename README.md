# RNA-Seq analyses

This is a basic introduction to using a computer cluster to carry out some analyses of RNA-Seq data.

# Getting set up

First we will ssh into pegasus

```bash
ssh USERNAME@pegasus.ccs.miami.edu
```

Next, we will download the files we need to carry out this exercise from Dropbox

```bash
wget --no-check-certificate https://www.dropbox.com/sh/enrzbdvvu3kluqe/AAAZJbGkPFob2tzsuFsDcToya?dl=1 -O taiwan_fq.zip
unzip taiwan_fq.zip
```

Next, we will build a STAR mapping index of our genome (here just one scaffold of the Pocillopora damicornis genome).
This will allow us to map reads from each of our samples

This script loads the STAR program onto the computer cluster, makes a directory to store the index, and runs STAR to build the index.
It looks like this:
```bash
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
```

OK. let's run it
```bash
bsub ./STARbuild.sh
```

# Mapping reads

This script is similar to the index building script, but instead of using STAR to build an index, loops through our samples and maps them to that index (I've omitted the BSUB header).
This will map each sample (each sample is a compressed file of raw reads i.e. a fq.gz file) to the genome to yield aligned read files (.bam files).
```bash
#!/bin/bash
#adapted from a script by Mike Connelly
#purpose: align trimmed RNAseq reads using STAR on the Pegasus bigmem queue

#BSUB -J staralign
#BSUB -U rsmasw

module load star

while read i; do
BASE=$(basename $i .fq.gz)
echo "Aligning $i" 
STAR \
--runMode alignReads \
--readFilesIn $i \
--readFilesCommand gunzip -c \
--genomeDir STARindex \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile pdam_1415_maker.gtf \
--outFileNamePrefix "$BASE"_unsorted
#lets me know file is done
echo "STAR alignment of $i complete"
#index the bam file
samtools view -b "$BASE"Aligned.out.sam > "$BASE"_unsorted.bam
samtools sort "$BASE"_unsorted.bam > "$BASE".bam
samtools index "$BASE".bam
rm "$BASE"Aligned.out.sam.sam
rm "$BASE"_unsorted.bam
done < samples.txt


```

let's map our reads!
```bash
#make a list of samples
ls -1 *.fq.gz > samples.txt
#submit job
bsub ./STARalign.sh
```

# Getting gene expression counts

Now we will count the number of reads that mapped to each gene in each sample

Let's look at the script
```bash
#!/bin/bash

#BSUB -J featurecounts
#BSUB -U rsmasw

module load subread

featureCounts -a pdam_1415_maker.gtf -o feature_counts.out -t exon -g gene_id bamlist.txt
```

Ok, let's run it!
```bash
ls -1 *.bam > bamlist.txt
bsub ./featureCounts.sh
```

# Retrieving data
We will now retrieve the expression data for analysis and visualization.
Open a new terminal without logging onto the cluster.
```bash
mkdir workshopAnalyses
cd workshopAnalyses
scp USERNAME@pegasus.ccs.miami.edu:feature_counts.out .
scp USERNAME@pegasus.ccs.miami.edu:WTsamples_all.txt .
```

# DESeq analysis

copy this into a new R script called basic_expression_analyses.R in Rstudio

```R
#basic_expression_analyses.R
###INPUT: metadata "WTsamples_all.txt" and gene expression counts "feature_counts.out"
###OUTPUT: normalized expression matrix "normCounts.csv", list of interesting genes "DEgenes.txt"

library('DESeq2')

#read in the metadata
meta<-read.delim('WTsamples_all.txt')

#make sure treatments are factors and not continuous variables for DESeq analysis
meta$colony<-factor(meta$colony)
meta$heat<-factor(meta$heat)
meta$antibiotics<-factor(meta$antibiotics)
meta$condition<-factor(meta$condition,levels=c('control','Heat','Antibiotics','Antibiotics.Heat','LPS','Antibiotics.Heat.LPS'))

#read in the counts
counts<-read.delim('feature_counts.out',skip=1)

#save the geneIDs so we can add them to our expression matrix
geneIDs<-as.character(counts[,1])

#remove extra gene information so we just have a matrix of gene expression counts
counts<-counts[,7:ncol(counts)]

#add gene names to expression matrix
rownames(counts)<-geneIDs

#make expression matrix sample names match our metadata
colnames(counts)<-meta$sample

#normalize the full expression matrix for library size and save for future plotting and analysis
cds<-DESeqDataSetFromMatrix(counts,meta,~colony+heat*antibiotics)
cds<-estimateSizeFactors(cds)
normCounts<-counts(cds,norm=T)
write.csv(normCounts,'normCounts.csv')

#subset data to just heat and antibiotic treatments for this analysis
counts<-counts[,meta$lps==0]
meta<-meta[meta$lps==0,]

#carry out DESeq analysis
cds<-DESeqDataSetFromMatrix(counts,meta,~colony+heat*antibiotics)
cds<-DESeq(cds)

#get DESeq results for heat treatment
res<-results(cds,name='heat_1_vs_0')

#get DESeq results for interaction between heat and antibiotics
res2<-results(cds,name='heat1.antibiotics1')

#save the names of the genes which responded to heat and had an interaction between antibiotics and heat response
sig<-res[which(res$padj<0.2&res2$padj<0.2),]
write.table(rownames(sig),quote=F,row.names=F,col.names=F,file='DEgenes.txt')
```

# Plotting

copy this into a new R script called plot_genes.R in Rstudio


```R

###INPUT: metadata "WTsamples_all.txt", normalized expression matrix "normCounts.csv", and list of interesting genes "DEgenes.txt"
###OUTPUT: plots of expression across treatments for each interesting gene

#read in metadata
meta<-read.delim('WTsamples_all.txt')

#order conditions to make looking at plots easier
meta$condition<-factor(meta$condition,levels=c('control','Heat','Antibiotics','Antibiotics.Heat','LPS','Antibiotics.Heat.LPS'))

#read in normalized expression matrix
normCounts<-as.matrix(read.csv('normCounts.csv',row.names=1))

#read in list of interesting genes from earlier analysis
genes<-read.delim('DEgenes.txt',header=F)[,1]

#plot each gene in a separate pdf
for(gene in genes){
	pdf(file=paste(gene,".pdf"))
	par(mar=c(8,4,2,1))
	stripchart(normCounts[gene,]~meta$condition,vertical=T,pch=1,method='j',las=2,ylab='Normalized counts',main=curr)
	dev.off()
}
```

# Annotation

we can also extract the interesting genes' protein sequences for further annotation. First we need to install samtools and get the protein sequences
```bash
#install local copy of samtools
wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
tar -xvjf samtools-1.8.tar.bz2
cd samtools-1.8
make
cd ..

#get protein sequences
scp USERNAME@pegasus.ccs.miami.edu:pdam_1415_maker.faa .
bash 
```

Then we can run the following script to extract the interesting genes' protein sequences.
Use a text editor to save the following in a file called get_DEgene_proteins.sh
```bash
#!/bin/bash
#USAGE: bash get_DEgene_proteins.sh DEgenes.txt proteins.faa output.faa

xargs samtools-1.8/samtools faidx $2 < $1 > $3
```

Let's run it (locally this time, so no bsub)!
```bash
bash ./get_DEgene_proteins.sh DEgenes.txt pdam_1415_maker.faa DEgenes.faa
```



