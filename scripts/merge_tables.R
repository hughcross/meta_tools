#!/usr/bin/env Rscript

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "haplotype_counts.tsv"
}

##
Personalized_Reader <- function(lambda){
  read.table(lambda, header = TRUE ,sep = '\t')}

ff <- list.files(path=args[1], full.names=TRUE)

myfilelist <- lapply(ff, Personalized_Reader)

names(myfilelist) <- list.files(path="/Users/hughcross/OneDrive/AnatAnalysis/Paua/june20_extractHaps/maps/logfiles/", full.names=FALSE)

bigtable <- myfilelist %>%
  reduce(left_join, by="Haplotype")

write.table(bigtable, file = args[2],sep="\t",row.names = FALSE, quote = FALSE)
