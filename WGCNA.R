setwd('~/workshopAnalyses')
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