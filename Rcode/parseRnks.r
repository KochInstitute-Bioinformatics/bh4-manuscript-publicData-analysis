library(openxlsx)
library(DESeq2)
library(tidyverse)

options(device=pdf)
setwd("Z:/charliew/da/gsea")

rnkFile <- read.xlsx("huMPG_neu_rnkStock.xlsx", colNames=TRUE, rowNames=TRUE)

names <- colnames(rnkFile)

for (i in names){
  d<-rnkFile %>% select(i)
  d$'#gene' <- row.names(d)
  d <- d[,c(2,1)]
  write.table(d, sep='\t',file=paste0(i,".rnk"),col.names=TRUE, quote=FALSE, row.names=FALSE)
}
