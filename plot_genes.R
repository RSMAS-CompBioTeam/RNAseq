###INPUT: metadata "WTsamples_all.txt", normalized expression matrix "normCounts.csv", and list of interesting genes "DEgenes.txt"
###OUTPUT: plots of expression across treatments for each interesting gene

setwd('~/workshopAnalyses')

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
