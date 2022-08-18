library(ggplot2)
library(dplyr)
library(ggrepel)

## load data
path="git/gitHub/RNA-isoforms-chek2/data/"
col.pal <- c(jessie="#DF7861", Hadley="#94B49F")

## Jessies data
jessie.bed <- read.delim(paste0(path,'chek2_jessieRun10.isoforms.bed'), header =F)
jessie.tsv <- read.delim(paste0(path,'chek2.jessieRun10.isoforms.tsv'))


# 5% threshold line
jess.thres <- quantile(jessie.tsv$count,probs = 0.05)
median.col <- 'black'
thres.col <- 'black'

p1 <- jessie.tsv %>% mutate(label=stringr::str_extract(ids,"ENST\\w+"), count=log2(count)) %>% 
  ggplot(aes(x=1,y=count, label=label)) + 
  geom_jitter(position = position_jitter(seed = 1), color=col.pal[1]) + 
  geom_hline(aes(yintercept = median(count)), color=median.col) + ylab("log2(counts)") +
  geom_hline(yintercept = log2(jess.thres), color=thres.col, linetype='dashed') + ylab("log2(counts)") +
  coord_flip()+ theme_bw() + geom_label_repel(position=position_jitter(seed=1),box.padding = 0.5,min.segment.length = unit(0, 'lines'))+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank())

p2 <- jessie.tsv %>% mutate(label=stringr::str_extract(ids,"ENST\\w+"), count=log2(count)) %>% 
  ggplot(aes(x=count, label=label)) + geom_histogram(bins=50,fill=col.pal[1], color= 'black') + xlab("Number of reads (log2)") + 
  ylab("Number of isoforms") +theme_bw() + 
  geom_vline(aes(xintercept = median(count)), color=median.col) + ylab("log2(counts)") +
  geom_vline(xintercept = log2(jess.thres), color=thres.col, linetype='dashed')
