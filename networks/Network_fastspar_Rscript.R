
library(dplyr)
library(magrittr)

path="pvalues.tsv"
pvals=read.csv2(path,header=FALSE,sep="\t")
rownames(pvals) <- pvals[,1]
colnames(pvals) <- pvals[1,]
pmat=pvals[2:nrow(pvals),2:ncol(pvals)]

path="median_correlation.tsv"
corvals=read.csv2(path,header=FALSE,sep="\t")
rownames(corvals) <- corvals[,1]
colnames(corvals) <- corvals[1,]
cormat=corvals[2:nrow(corvals),2:ncol(corvals)]

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
network=flattenCorrMatrix(cormat, pmat)
as.data.frame(network)
network= network %>% filter(cor > 0.6 | cor > -0.6,
                            p < 0.05)
write.table(network,file="netcyt.txt",sep="\t", quote=FALSE)

