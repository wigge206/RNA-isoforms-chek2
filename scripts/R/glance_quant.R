## Quick summary of reads mapped to genes
## Load packages
library(ggplot2)
library(dplyr)
library(cowplot)

## global paramaters
col.pal <- c(jessie="#DF7861", Hadley="#94B49F")

## load data
setwd("/WORKSPACE/George/CHEK2_RNAisoforms/quantify")
jessie.dat <- read.delim("jessieRun10_countMatrix.tsv") ## Jessie's data 
hadley.dat <- read.delim("hadley_barcode91_95_countMatrix.tsv") ## Jessie's data

## extract genes - recode to hgnc symbols - remove reads not mapped to genes
jessie.dat$geneID <- stringr::str_extract(jessie.dat$ids, "ENSG\\w+")
hadley.dat$geneID <- stringr::str_extract(hadley.dat$ids, "ENSG\\w+")
jessie.dat <- jessie.dat %>% filter(!is.na(geneID)) %>% mutate(geneID = recode(geneID, 
'ENSG00000183765' ='CHEK2',
'ENSG00000136492' = 'BRIP1',
'ENSG00000108384' = 'RAD51C',
'ENSG00000185379' = 'RAD51D',
.default = NA_character_  ## all other genes as NA
))

hadley.dat <- hadley.dat %>% filter(!is.na(geneID)) %>% mutate(geneID = recode(geneID, 
'ENSG00000183765' ='CHEK2',
'ENSG00000039068' = 'CDH1',
'ENSG00000083093' = 'PALB2',
'ENSG00000171862' = 'PTEN',
.default = NA_character_  ## all other genes as NA
))

## Plot reads mapped to each gene
p1 <- jessie.dat %>% group_by(geneID) %>% summarise(count=sum(sample1_condition1_batch1)) %>% 
ggplot(aes(reorder(geneID, -count), count)) + geom_col(fill=col.pal[1]) + theme_bw()+ xlab(NULL) + ylab("Number of reads")

p2 <- hadley.dat %>% group_by(geneID) %>% summarise(across(contains("sample"),sum)) %>% mutate(count= rowSums(select(., contains('sample')))) %>%
ggplot(aes(reorder(geneID, -count), count)) + geom_col(fill=col.pal[2]) + theme_bw()+ xlab(NULL) + ylab("Number of reads")
title <- ggdraw() +draw_label("Number of reads mapped to genes",
    fontface = 'bold', x = 0, hjust = 0) +theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(title, plot_grid(p1,p2, labels=LETTERS[1:2]), ncol=1, rel_heights=c(0.1,1)) 
ggsave(file="reads_Summary_geneLevel.png", width=7, height=3)

## plot chek2 RNA isoforms
jessie.dat.chek2 <- jessie.dat %>% filter(geneID == "CHEK2") %>% arrange(desc(sample1_condition1_batch1)) %>% 
    mutate(id2= ifelse(sample1_condition1_batch1 >= 2273, ids, 'other')) %>% group_by(id2) %>% summarise(count= sum(sample1_condition1_batch1))
p3<- jessie.dat.chek2  %>% ggplot(aes(reorder(id2,-count), count)) +theme_bw()+ 
geom_col(fill=col.pal[1]) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+xlab(NULL) + ylab("Number of reads")

hadley.dat.chek2 <- hadley.dat %>% filter(geneID == "CHEK2") %>% tidyr::pivot_longer(cols=2:6)
chek2.means <- hadley.dat %>% filter(geneID == "CHEK2") 
chek2.means$avg <- rowMeans(chek2.means[,2:6])
top10 <- chek2.means %>% arrange(desc(avg)) %>% pull(ids) %>% head(10)
p4 <- hadley.dat.chek2 %>% mutate(id2 = ifelse(ids %in% top10, ids, 'other')) %>% group_by(id2,name) %>% summarise(count= sum(value)) %>%
mutate(name = recode(name,
'sample1_conditionA_batch1' = 'Sample 1',
'sample2_conditionA_batch1' = 'Sample 2',
'sample3_conditionA_batch1' = 'Sample 3',
'sample4_conditionA_batch1' = 'Sample 4',
'sample5_conditionA_batch1' = 'Sample 5',
)) %>% filter(name != 'Sample 5') %>% ## has no reads
ggplot(aes(reorder(id2,-count), count)) + geom_col(fill=col.pal[2]) + facet_wrap(.~ name, scales='free')  +theme_bw() +
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +xlab(NULL) + ylab("Number of reads")

## multiplot-title
title <- ggdraw() +draw_label("Number of reads mapped to CHEK2 isoforms",
    fontface = 'bold', x = 0, hjust = 0) +theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(title,plot_grid(p3,p4, nrow=1, labels=LETTERS[1:2], rel_widths=c(.4,.6)), ncol=1, rel_heights=c(.1,1))
ggsave(file="reads_Summary_isoFormlevel.png", width=7, height=3)