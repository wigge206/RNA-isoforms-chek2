#############
# Description: I have manually annotated a bunch of isoform for hadley - I will copy these to Jessie
############
library(stringr)
library(dplyr)


jessie <- read.delim("git/gitHub/RNA-isoforms-chek2/data/chek2_jessie_annotated.txt")
hadley <- read.delim("git/gitHub/RNA-isoforms-chek2/data/chek2_hadley_annotated.txt")


j.df <- data.frame(ids = jessie$id, 
                   exon2match = jessie$exons-1, 
                   exonSize2match = gsub("^\\w+,","",gsub("\\w+,\\w+,$","",jessie$exonSizes)),
                   exonStart2match = gsub("^\\w+,","",gsub("\\w+,$","",jessie$exonStarts)),
                   counts= jessie$Counts
                   )
## Need to adjust exonStart2match. Jessies reads map 153 bp upstream of hadleys
j.df$exonStart2match <-unlist(lapply(strsplit(j.df$exonStart2match, ","), function(i) paste(as.numeric(i)- 153, collapse=",")))

h.df <- data.frame(id = hadley$id,
                   exon2match = hadley$exons,
                   exonSize2match = gsub("^\\w+,","",gsub("\\w+,$","",hadley$exonSizes)),
                   exonStart2match =gsub("^\\w+,","", gsub(",$","",hadley$exonStarts)),
                   Events = hadley$Events,
                   r.ENST00000404276.6 = hadley$r...ENST00000404276.6
                   )



test<- join(j.df, h.df, by=c('exon2match','exonSize2match','exonStart2match')) 

write.csv(test, "tmp.csv")
