
library(SummarizedExperiment)
library(ggplot2)
library(DESeq2)
library(HandyPack)
library(tictoc)
library(stringr)
library(pheatmap)
library(BioHandy)


rm(list=ls())
graphics.off()
## ###################################################
## ###################################################
addReplicate = function(se)
{
    colData = colData(se)

    a = str_split(colData$SampleName,'_')
    colData$replicate = unlist(lapply(a,function(x) return(x[1])))
    colData(se) = colData

    return(se)
}

## ###################################################
rowToDF = function(se,row)
{
    cd = colData(se)
    df = data.frame(Y=assay(se)[row,],
                    replicate=cd$replicate,
                    cells=cd$cells,
                    treatment=cd$treatment)

    return(df)
}

## ###################################################
regressOutReplicate = function(se)
{
    M = matrix(0,nrow=nrow(se),ncol=ncol(se))
    rownames(M) = rownames(se)
    colnames(M) = colnames(se)
    
    for(i in 1:nrow(se))
    {
        df = rowToDF(se,i)
        m = lm(Y~replicate,df)
        M[i,] = df$Y - predict(m)
    }
    return(M)
}

## ###################################################

## Load the SE:
Tic('loading se')
load('grandSE/grandSE.RData')
grandSE = addReplicate(grandSE)
toc()

## Omit unexpressed genes:
idx = rowSums(assay(grandSE)) > 0
grandSE = grandSE[idx,]
grandSE = addReplicate(grandSE)

## Get gene names:
Tic('genes')
genes = entrezToGeneSymbol(rownames(grandSE),'human')
toc()
genes = genes[!is.na(genes)]
genes = genes[genes != 'NA']

for(i in 1:nrow(grandSE))
    if(rownames(grandSE)[i] %in% names(genes))
        rownames(grandSE)[i] = genes[rownames(grandSE)[i]]


dds = DESeqDataSet(grandSE,design=~ replicate + treatment)
vStab = varianceStabilizingTransformation(dds)


Tic('regressing')
M = regressOutReplicate(vStab)
toc()

stopifnot(identical(colnames(M),rownames(colData(grandSE))))
colnames(M) = colData(grandSE)$SampleName

idx = str_detect(colnames(M),'Diff_DMSO') |
    str_detect(colnames(M),'Diff_CH')

diffs = M[,idx]
diffs = diffs[,c(3,5,1,4,6,2)]

Tic('pheatmap')
jpeg(filename='replicateHeatmap.jpg',
height=20,width=12,units='in',res=100)
pheatmap(diffs,
         cluster_rows=TRUE,
         cluster_cols=FALSE,
         treeheight_row=0,
         treeheight_col=0,
         show_rownames=FALSE)
dev.off()
toc()


stdev = numeric(nrow(diffs))
for(i in 1:nrow(diffs))
    stdev[i] = sd(diffs[i,])

anotherM = diffs
anotherM = anotherM[order(-stdev),]
anotherM = anotherM[1:100,]

jpeg(filename='replicateHeatmapTop100.jpg',
height=15,width=12,units='in',res=100)
pheatmap(anotherM,
         cluster_rows=TRUE,
         cluster_cols=FALSE,
         treeheight_row=0,
         treeheight_col=0)
dev.off()

