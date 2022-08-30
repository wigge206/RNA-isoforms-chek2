library(dplyr)


path="git/gitHub/RNA-isoforms-chek2/data/"
jessie.dat <- read.delim(paste0(path,"chek2.jessieRun10.isoforms.tsv"))  
hadley.dat <- read.delim(paste0(path,"chek2.hadley91_95isoforms.tsv")) 
jessie.bed <- read.delim(paste0(path,'chek2_jessieRun10.isoforms.bed'), header =F)
hadley.bed <- read.delim(paste0(path,'chek2_hadley91_95.isoforms.bed'), header =F)


## top 10 files to look at
hadley.dat$TotalSum <- rowSums(hadley.dat[,2:6])
ids_had <- hadley.dat %>% slice_max(TotalSum, n=10) %>% pull(ids)
ids_jess<-jessie.dat %>% slice_max(count, n=10) %>% pull(ids)

sub_had <- hadley.bed[match(ids_had,hadley.bed$V4),]
sub_had$V4 <-  gsub("_ENSG00000183765.23" ,"", sub_had$V4)
j = 1
for(i in 1:nrow(sub_had)){
  if(!grepl("ENST", sub_had$V4[i])){
    sub_had$V4[i] <- paste0("FLAIR.tx.",j)
    j=j+1
  }
}

sub_jess <- jessie.bed[match(ids_jess,jessie.bed$V4),]
sub_jess$V4 <-  gsub("_ENSG00000183765.23" ,"", sub_jess$V4)
j = 1
for(i in 1:nrow(sub_jess)){
  if(!grepl("ENST", sub_jess$V4[i])){
    sub_jess$V4[i] <- paste0("FLAIR.tx.",j)
    j=j+1
  }
}

write.table(file="FilesForVisualisationONLY/top10_hadley.bed", sub_had, row.names=F, col.names=F, quote=F)
write.table(file="FilesForVisualisationONLY/top10_Jessie.bed", sub_jess, row.names=F, col.names=F, quote=F)

