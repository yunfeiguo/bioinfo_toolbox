#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
table <- args[2]
load(file)
write.table(eval(parse(text=table)),file=paste(file,".csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
