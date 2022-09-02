library(dplyr)
library(ggplot2)
library(openxlsx)
## global parameters
col.pal <- c(jessie="#DF7861", Hadley="#94B49F")

path="git/gitHub/RNA-isoforms-chek2/data/"

jessie.bed <- read.delim(paste0(path,'chek2_jessieRun10.isoforms.bed'), header =F)
hadley.bed <- read.delim(paste0(path,'chek2_hadley91_95.isoforms.bed'), header =F)
gencode.bed <- read.delim("chek2.gencode.41.annotation.bed", header=F)
# MANE transcript is ENST00000404276.6

bind_rows(jessie.bed,hadley.bed, .id = "id" )%>% 
  mutate(id=recode_factor(id, '1'="Jessie's",'2'="Hadley's")) %>% 
  ggplot(aes(V10))+ geom_histogram(aes(fill=id), binwidth=1,color='white')+  theme_bw() +
  scale_fill_manual(values=as.character(col.pal))+ ylab("Count") + xlab(NULL)+
  facet_wrap(.~id) +theme(legend.position = 'none') +
  scale_x_continuous(breaks = seq(4,18,2), lim = c(3,19)) + ggtitle("Exon per isoform frequency")
ggsave(file="git/git-site/RNA-isoforms-chek2/assets/images/chek2ExonFreq.png", width = 7, height=3)


## load in counts to order
jessie.tsv <- read.delim(paste0(path,'chek2.jessieRun10.isoforms.tsv'))
hadley.tsv <- read.delim(paste0(path,'chek2.hadley91_95isoforms.tsv'))

jessie.tsv <- jessie.tsv %>% arrange(desc(count))
id.order <- jessie.tsv %>% arrange(desc(count)) %>% pull(ids)
jessie.bed <- jessie.bed[match(id.order, jessie.bed$V4),]

hadley.tsv$count = rowSums(hadley.tsv[,2:6])
hadley.tsv <- hadley.tsv %>% arrange(desc(count))
id.order <- hadley.tsv %>% arrange(desc(count)) %>% pull(ids)
hadley.bed <- hadley.bed[match(id.order, hadley.bed$V4),]

## Isoform renaming - the long readName_geneName ids are annoying
## Rename as FLAR.tx.# for unknown, for known strip ENSG id

## generate and ID lookup in case I need to remap.
id_lookup.jessie <- data.frame(orginial =jessie.bed$V4, gw_annotation = jessie.bed$V4)
id_lookup.hadley <- data.frame(orginial =hadley.bed$V4, gw_annotation = hadley.bed$V4)
id_lookup.jessie$gw_annotation<-gsub("_ENSG\\w+.\\w+$","",id_lookup.jessie$gw_annotation)
id_lookup.hadley$gw_annotation<-gsub("_ENSG\\w+.\\w+$","",id_lookup.hadley$gw_annotation)

j = 1
for(i in 1:nrow(id_lookup.jessie)){
  if(!grepl("ENST", id_lookup.jessie$gw_annotation[i])){
    id_lookup.jessie$gw_annotation[i] <- paste0("FLAIR.tx.",j)
    j=j+1
  }
}
j = 1
for(i in 1:nrow(id_lookup.hadley)){
  if(!grepl("ENST", id_lookup.hadley$gw_annotation[i])){
    id_lookup.hadley$gw_annotation[i] <- paste0("FLAIR.tx.",j)
    j=j+1
  }
}


####### Save ID look up
wb <- createWorkbook()
addWorksheet(wb, sheetName = "Jessie_id")
addWorksheet(wb, sheetName = "hadley_id")

writeData(wb,"Jessie_id", id_lookup.jessie)
writeData(wb,"hadley_id", id_lookup.hadley)

saveWorkbook(wb, file = "git/gitHub/RNA-isoforms-chek2/data/isoform_idlookup.xlsx", overwrite = TRUE)
#######

jessie.bed$V4<-id_lookup.jessie$gw_annotation
hadley.bed$V4<-id_lookup.hadley$gw_annotation

jessie.bed$counts <- jessie.tsv$count
hadley.bed <- cbind(hadley.bed, hadley.tsv[,2:7])


write.table(jessie.bed, "git/gitHub/RNA-isoforms-chek2/data/chek2_jessie_annotated.txt", row.names =F, col.names=F, sep="\t")
write.table(hadley.bed, "git/gitHub/RNA-isoforms-chek2/data/chek2_hadley_annotated.txt", row.names =F, col.names=F, sep="\t")
