library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

#### load data
path="git/gitHub/RNA-isoforms-chek2/data/"
col.pal <- c(jessie="#DF7861", Hadley="#94B49F")

## Jessies data
jessie.bed <- read.delim(paste0(path,'chek2_jessieRun10.isoforms.bed'), header =F)
jessie.tsv <- read.delim(paste0(path,'chek2.jessieRun10.isoforms.tsv'))
## Hadley's data
hadley.bed <- read.delim(paste0(path,'chek2_hadley91_95.isoforms.bed'), header =F)
hadley.tsv <- read.delim(paste0(path,'chek2.hadley91_95isoforms.tsv'))

hadley.tsv$total <- rowSums(hadley.tsv[2:6])


### 5% threshold line
jess.thres <- quantile(jessie.tsv$count,probs = 0.9)
hadley.thres <- quantile(hadley.tsv$total,probs = 0.1)

### Color palate
median.col <- '#6189df'
thres.col <- '#6189df'

p1 <- jessie.tsv %>% mutate(label=stringr::str_extract(ids,"ENST\\w+"), count=log2(count)) %>% 
  ggplot(aes(x=1,y=count, label=label)) + 
  geom_jitter(position = position_jitter(seed = 1), color=col.pal[1]) + ylim(0,20)+ 
  geom_hline(aes(yintercept = median(count)), color=median.col,size=1) + ylab("log2(counts)") +
  geom_hline(yintercept = log2(jess.thres), color=thres.col, linetype='dashed', size=1) +
  coord_flip()+ theme_bw() + 
  geom_label_repel(position=position_jitter(seed=1),box.padding = 0.5,min.segment.length = unit(0, 'lines'), size=3)+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank())

p2 <- jessie.tsv %>% mutate(label=stringr::str_extract(ids,"ENST\\w+"), count=log2(count)) %>% 
  ggplot(aes(x=count, label=label)) + geom_histogram(bins=50,fill=col.pal[1], color= 'black') + xlab("Number of reads (log2)") + 
  ylab("Number of isoforms") +theme_bw() + xlim(0,20)+
  geom_vline(aes(xintercept = median(count)), color=median.col, size=1) + ylab("log2(counts)") +
  geom_vline(xintercept = log2(jess.thres), color=thres.col, linetype='dashed', size=1)

title <- ggdraw() +draw_label("Jessie's Data: Distribution of reads",
                              fontface = 'bold', x = 0, hjust = 0) +theme(plot.margin = margin(0, 0, 0, 7))

top_row <- plot_grid(p1,p2)
top_row <- plot_grid(title, top_row,ncol=1, rel_heights=c(0.1,1))

p3 <- hadley.tsv %>% mutate(label=stringr::str_extract(ids,"ENST\\w+"), count=log2(total)) %>%ggplot(aes(x=1,y=count, label=label)) + 
  geom_jitter(position = position_jitter(seed = 1), color=col.pal[2]) + ylim(0,20)+ 
  geom_hline(aes(yintercept = median(count)), color=median.col,size=1) + ylab("log2(counts)") +
  geom_hline(yintercept = log2(hadley.thres), color=thres.col, linetype='dashed', size=1) +
  coord_flip()+ theme_bw() + 
  geom_label_repel(position=position_jitter(seed=1),box.padding = 0.5,min.segment.length = unit(0, 'lines'),size=3)+
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank())

p4 <- hadley.tsv %>% mutate(label=stringr::str_extract(ids,"ENST\\w+"), count=log2(total)) %>% 
  ggplot(aes(x=count, label=label)) + geom_histogram(bins=50,fill=col.pal[2], color= 'black') + xlab("Number of reads (log2)") + 
  ylab("Number of isoforms") +theme_bw() + xlim(0,20)+
  geom_vline(aes(xintercept = median(count)), color=median.col, size=1) + ylab("log2(counts)") +
  geom_vline(xintercept = log2(hadley.thres), color=thres.col, linetype='dashed', size=1)

title <- ggdraw() +draw_label("Hadley's Data: Distribution of reads",
                              fontface = 'bold', x = 0, hjust = 0) +theme(plot.margin = margin(0, 0, 0, 7))

bottom_row <- plot_grid(p3,p4)
bottom_row <- plot_grid(title, bottom_row,ncol=1, rel_heights=c(0.1,1))


plot_grid(top_row,bottom_row,ncol=1)
ggsave("git/git-site/RNA-isoforms-chek2/assets/images/Read_isoformsDistrubtion.png", width=7, height=5)
