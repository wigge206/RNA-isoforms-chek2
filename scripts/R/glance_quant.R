## Quick summary of reads mapped to genes
## Load packages
library(ggplot2)
library(dplyr)
library(cowplot)

## load data
setwd("/WORKSPACE/George/CHEK2_RNAisoforms/quantify")
dat <- read.delim("jessieRun10_countMatrix.tsv") ## Jessie's data 
dat$geneID <- stringr::str_extract(dat$ids, "ENSG\\w+")

p1 <- dat %>% filter(!is.na(geneID)) %>% group_by(geneID) %>% summarise(count=sum(sample1_condition1_batch1)) %>% 
mutate(gene = ifelse(count > 10000, geneID, "other"), 
gene = recode(gene, 
'ENSG00000183765' ='CHEK2',
'ENSG00000136492' = 'BRIP1',
'ENSG00000108384' = 'RAD51C',
'ENSG00000185379' = 'RAD51D'
)) %>% ggplot(aes(reorder(gene,-count), count)) + geom_col(fill="#383d57") +theme_bw()+ xlab(NULL) + ylab("Number of reads") +
ggtitle("Reads per gene")

p2<-dat  %>% rename(count= 'sample1_condition1_batch1') %>%
mutate(gene = recode(geneID, 
'ENSG00000183765' ='CHEK2',
'ENSG00000136492' = 'BRIP1',
'ENSG00000108384' = 'RAD51C',
'ENSG00000185379' = 'RAD51D',
.default = NA_character_
))%>% filter(!is.na(gene)) %>% group_by(gene) %>% mutate(labels=ifelse(count/sum(count) < 0.01, 'other', ids)) %>%
 ggplot(aes(x=reorder(labels, -count), count)) + theme_bw()+
geom_col(color="#383d57",fill="#383d57") + facet_wrap(.~gene, scales='free') + 
theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+ xlab(NULL) + ylab("Number of reads") +ggtitle("Reads per isoform")

png(file="reads_Summary_geneLevel.png", width=10.9, height=3.56,units='in', res=300)
plot_grid(p1,p2, nrow=1)
dev.off()


## Hadley's data
dat <- read.delim("hadley_barcode91_95_countMatrix.tsv") ## Jessie's data 
dat$geneID <- stringr::str_extract(dat$ids, "ENSG\\w+")
dat$avg <- rowMeans(dat[,2:6])

dat<-tidyr::pivot_longer(dat, cols=2:6)

x[rev(order(x$`sum(value)`)),]




gene = recode(gene, 
'ENSG00000183765' ='CHEK2',
'ENSG00000039068' = 'CDH1',
'ENSG00000083093' = 'PALB2',
'ENSG00000171862' = 'PTEN'
))



"#57384c"