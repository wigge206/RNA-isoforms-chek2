library(dplyr)
library(ggplot2)

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
