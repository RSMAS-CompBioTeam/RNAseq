library('DESeq2')

meta<-read.delim('WTsamples_all.txt')
counts<-read.delim('~/Dropbox/RSMAS-CompBio-athon-2018/feature_counts.out',skip=1)
geneLengths<-counts$Length
geneIDs<-as.character(counts[,1])
counts<-counts[,7:ncol(counts)]
rownames(counts)<-geneIDs
colnames(counts)<-meta$sample
meta$condition<-factor(meta$condition,levels=c('control','Heat','LPS','Antibiotics','Antibiotics.Heat','Antibiotics.Heat.LPS'))

cds<-DESeqDataSetFromMatrix(counts,meta,~colony+condition)
cds<-DESeq(cds,test='LRT',reduced=~colony)
res<-results(cds)
hist(res$pvalue)
sig<-res[which(res$padj<0.01),]
write.table(rownames(sig),quote=F,row.names=F,col.names=F,file='DEgenes.txt')

normCounts<-counts(cds,norm=T)
write.csv(normCounts,'normCounts.csv')
#stripchart(normCounts[rownames(sig)[35],]~meta$condition,vertical=T,method='j',pch=1,las=2,cex.axis=0.5)