# RNA-Seq analyses

This is a basic introduction to using a computer cluster to carry out some analyses of RNA-Seq data.

# Getting set up

First we will ssh into pegasus

```bash
ssh USERNAME@pegasus.ccs.miami.edu
bsub -Is -P ccsfellows bash
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
#BSUB -P ccsfellows
#BSUB -o %J.out

module load star

mkdir STARindex

STAR \
--runMode genomeGenerate \
--genomeDir STARindex \
--genomeFastaFiles pdam_genome_1415.fasta \
--sjdbGTFfile pdam_1415_maker.gtf
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
#BSUB -P ccsfellows
#BSUB -o %J.out

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
--outFileNamePrefix "$BASE"
#lets me know file is done
echo "STAR alignment of $i complete"
#index the bam file
samtools view -b "$BASE"Aligned.out.sam > "$BASE"_unsorted.bam
samtools sort "$BASE"_unsorted.bam > "$BASE".bam
samtools index "$BASE".bam
rm "$BASE"Aligned.out.sam
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

#BSUB -J staralign
#BSUB -U rsmasw
#BSUB -P ccsfellows
#BSUB -o %J.out


module load subread

featureCounts -a pdam_1415_maker.gtf -o feature_counts.out -t exon -g gene_id Wt*.bam
```

Ok, let's run it!
```bash
bsub ./featureCounts.sh
```

# Retrieving data
We will now retrieve the expression data for analysis and visualization.
Open a new terminal without logging onto the cluster. We will copy the gene expression counts data (feature_counts.out) and the sample metadata (WTsamples_all.txt) into the CompBio-RNAseq Github repository we created earlier.

```bash
# Navigate to the CompBio-RNAseq repository
cd ~/github/CompBio-RNAseq

# Transfer the files from Pegasus into the data subdirectory
scp USERNAME@pegasus.ccs.miami.edu:~/feature_counts.out .
scp USERNAME@pegasus.ccs.miami.edu:~/WTsamples_all.txt .
```

# Open RStudio
The differential expression analysis will take place using R, so open RStudio. We will create a new Rproject associated with our git repository, which will help us stay organized and interact with Github.com from inside RStudio.

Go to File > New Project  
Click Existing Directory  
Choose your CompBio-RNAseq folder  

Now we will see how to add, commit, and push within the RStudio environment...

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
	stripchart(normCounts[gene,]~meta$condition,vertical=T,pch=1,method='j',las=2,ylab='Normalized counts',main=gene)
	dev.off()
}
```

# Annotation

we can also extract the interesting genes' protein sequences for further annotation.
first we transfer our list of genes to the cluster
```bash
scp DEgenes.txt USERNAME@pegasus.ccs.miami.edu:~/ 
```

Then we can run the following command to extract the interesting genes' protein sequences.
```bash
module load samtools

xargs samtools faidx pdam_1415_maker.faa < DEgenes.txt > DEgenes.faa
```

Now that we have the sequences of our differentially expressed genes, we can compare them to known sequences in a database to get information about these genes. 
Let's download the UniProt SwissProt database and BLAST our query sequences against it.

```bash
# Make a new directory in your nethome for the SwissProt database
mkdir ~/swissprot && cd ~/swissprot

# Download and unzip the SwissProt database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Make a BLAST database out of the SwissProt fasta file
module load blast
makeblastdb -in uniprot_sprot.fasta -dbtype prot
```

Now we have a BLAST database ready to compare our query sequences to. We will use the command `blastp` since we are comparing protein sequences against a protein database.

```bash
# Change back to home directory (where the DEgenes.faa file should be)
cd ~

# Run blastp with options
bsub -J blast -P ccsfellows \
blastp -query DEgenes.faa \
-db ~/swissprot/uniprot_sprot.fasta \
-out sprot_blastout.txt \
-evalue 1e-10
```

Now we can inspect the blast results (`less sprot_blastout.txt`) to see what proteins our differentially expressed genes are similar to. The default blastp output format shows all of the hits and the alignments between the query and subject sequences. If we want a more condensed output, there are many different formats we can choose from (see `blastp -help`). Let's try `-outfmt 6` to get a tabular output and add the argument `-num_alignments 1` to give us just the top hit fot each query.

```bash
# Run blastp with options
bsub -J blast -P ccsfellows \
blastp -query DEgenes.faa \
-db ~/swissprot/uniprot_sprot.fasta \
-out sprot_blastout_tabular.txt \
-evalue 1e-10 \
-outfmt 6 \
-num_alignments 1
```

Now we can inspect this tabular output (`less sprot_blastout_tabular.txt`). This format is well-suited to import into R for downstream analysis.

Search the UniProt website using the UniProt ID's of the top blast hits (e.g., D9IQ16) to get more information about these proteins.

# WGCNA
if you have a really large number of genes (which you often do in an RNA seq dataset) you might cluster their expression profiles
this can make analysis and interpretation easier

```bash
scp USERNAME@pegasus.ccs.miami.edu:~/meta_WGCNA.csv data/
scp USERNAME@pegasus.ccs.miami.edu:~/norm_counts_WGCNA.csv data/
```


```R
library(WGCNA)

counts<-read.csv('norm_counts_WGCNA.csv',row.names=1)
meta<-read.csv('meta_WGCNA.csv',row.names=1)
counts<-counts[,match(rownames(meta),colnames(counts))]
datExpr<-t(counts)

#find correlation power R^N that satisfies scale free critereon (SFT.R.sq>0.8)
sft<-pickSoftThreshold(datExpr,verbose=5)

#cluster genes into coexpressed modules, this may take a few minutes, or if you want it to run faster, uncomment the following line to filter for expression:
#datExpr<-datExpr[,apply(datExpr,2,mean)>150]
net<-blockwiseModules(datExpr,power=sft$powerEstimate,verbose=5,numericLabels=T)

#get "eigengenes" which summarize expression of all genes in a coexpressed module
MEs<-net$MEs

#see how tightly correlated individual genes are with their eigengene
membership<-cor(datExpr,MEs)

#gather what we know about the genes
geneInfo<-data.frame(row.names=colnames(datExpr),module=net$colors,membership=NA)
for(gene in rownames(geneInfo)){
	currmod<-geneInfo[gene,'module']
	geneInfo[gene,'membership']<-membership[gene,paste('ME',currmod,sep='')]
}
geneInfo<-geneInfo[order(geneInfo$module,geneInfo$membership,decreasing=T),]



#plot eigengenes along with highly correlated individual genes
pdf(file='module1.pdf',width=4,height=4)

#initialize plot with eigengene -- we scale the expression value to mean=0 sd=1 so we can compare profiles of genes with different absolute expression levels
MEagg<-aggregate(scale(MEs[,'ME1'])~Treatment+Time,FUN=mean,data=meta)
plot(MEagg[MEagg$Treatment=='H',2],(MEagg[MEagg$Treatment=='H',3]),type='l',lwd=3,col='red',ylim=c(-1,3),xlab='Time (minutes)',ylab='Scaled expression')

#plot each individual gene with high membership on top of the eigengene
for(gene in rownames(geneInfo)[geneInfo$module==1&geneInfo$membership>0.9]){
	agg<-aggregate(scale(datExpr[,gene])~Treatment+Time,FUN=mean,data=meta)
	lines(agg[agg$Treatment=='H',2],(agg[agg$Treatment=='H',3]),type='l',lwd=0.2,col='darkred')
	lines(agg[agg$Treatment=='C',2],(agg[agg$Treatment=='C',3]),type='l',lwd=0.2,col='darkblue')	
}

#replot eigengenes on top of the individual gene expression profiles
lines(MEagg[MEagg$Treatment=='H',2],(MEagg[MEagg$Treatment=='H',3]),type='l',lwd=3,col='red')
lines(MEagg[MEagg$Treatment=='C',2],(MEagg[MEagg$Treatment=='C',3]),type='l',lwd=3,col='blue')
dev.off()
```
