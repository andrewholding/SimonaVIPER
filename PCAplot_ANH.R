dds <- makeExampleDESeqDataSet(betaSD=1)
rld <- rlog(dds)
plotPCA(rld)

# also possible to perform custom transformation:
dds <- estimateSizeFactors(dds)
# shifted log of normalized counts
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1),
                           colData=colData(dds))
# the call to DESeqTransform() is needed to
# trigger our plotPCA method.
plotPCA( DESeqTransform( se ) )



#Run after the markdown. 

ddsMCF7<-DESeqDataSetFromMatrix(countData = MatrixMCF7, colData = coldata, 
                                design= ~Gene)

ddsMCF7 <- estimateSizeFactors(ddsMCF7)

genesMCF= read.csv("genesMCF7.txt", header= FALSE, sep=',')
splitted<- str_split_fixed(genesMCF$V2, "_",2)
splitted2<- str_split_fixed(splitted[,2], "_",2)
rownames(splitted2)<-genesMCF$V1

coldata_rep<-cbind(coldata,splitted2[rownames(coldata),2])
colnames(coldata_rep)[2]<-"Rep"

seMCF7<-SummarizedExperiment(log2(counts(ddsMCF7, normalized=TRUE)+1),colData =coldata_rep )

plotPCA( DESeqTransform(seMCF7),intgroup="Rep" )

