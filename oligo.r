#!/usr/bin/env Rscript

# loads package limma
library(oligo)

# substitute platform_name with the name of the platform used
platform_name="GPL57"

# read names of cel files
list=list.files(full.names=TRUE)

# read cel files
data=read.celfiles(list)

# background correction, normalization and probe summarization using RMA method
data.rma=rma(data)
data.matrix=exprs(data.rma)

# write the normalized expression in a text file
write.table(data.matrix,file=paste("data",platform_name,".txt",sep=""))
