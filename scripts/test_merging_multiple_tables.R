

setwd("/Users/hughcross/OneDrive/AnatAnalysis/Paua/june20_extractHaps/maps/logfiles/")

library(dplyr)
library(tidyverse)
## read in list of samples
samplelist <- scan('/Users/hughcross/OneDrive/AnatAnalysis/Paua/june20_extractHaps/maps/samplelist', what="", sep="\n")

samplelist
## loop through samples to read in tables
## https://stackoverflow.com/questions/22689301/how-use-read-table-in-a-for-loop

ff <- list.files(path="/Users/hughcross/OneDrive/AnatAnalysis/Paua/june20_extractHaps/maps/logfiles/", full.names=TRUE)
ff

## add a personalised reader
Personalized_Reader <- function(lambda){
  read.table(lambda, header = TRUE ,sep = '\t')}

myfilelist <- lapply(ff, Personalized_Reader)
names(myfilelist) <- list.files(path="/Users/hughcross/OneDrive/AnatAnalysis/Paua/june20_extractHaps/maps/logfiles/", full.names=FALSE)
myfilelist$`Sample-10_2_1.log`
samplelist[[1]]

## join 
typeof(myfilelist)

test <- left_join(myfilelist$`Sample-10_2_1.log`,myfilelist$`Sample-10_2_4.log`, by="Haplotype")
test

## try purr
## https://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list

bigtable <- myfilelist %>%
  reduce(left_join, by="Haplotype")

View(bigtable)

# output 
write.table(bigtable, file = "trial_haplotype_counts.txt",sep="\t",row.names = FALSE, quote = FALSE)










