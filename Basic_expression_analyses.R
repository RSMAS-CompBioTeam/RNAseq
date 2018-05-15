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