#!/usr/bin/env Rscript
library(LaplacesDemon)
library("optparse")
 
option_list = list(
  make_option(c("-f", "--logfile"), type="character", default=NULL, 
              help="name of the revbayes log file to use", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="name of the outdir to write ESS.txt to", metavar="character"),
  make_option(c("-b", "--burnin"), type="numeric", default=1000, 
              help="burnin iterations", metavar="numeric")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
df=read.table(opt$logfile, header=T, sep="\t")
df=df[-1:-opt$burnin,]

file=file.path(opt$outdir,"ESS.txt")
write("PARAMETER_NAME\tESS", file=file, append=T)
for (i in 2:ncol(df)) write(paste(colnames(df)[i], "\t", round(ESS(df[,i]),0), sep=""), file=file, append=T)
