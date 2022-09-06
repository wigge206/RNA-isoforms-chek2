library(dplyr)
library(openxlsx)
jessie.dat <- read.delim("git/gitHub/RNA-isoforms-chek2/data/chek2_jessie_annotated.txt")
hadley.dat <- read.delim("git/gitHub/RNA-isoforms-chek2/data/chek2_hadley_annotated.txt")



## Tidy up
total.counts <- sum(jessie.dat$Counts)
j.annotated<-jessie.dat %>% filter(Events !="", is.na(not.check)) %>% group_by(Events,r.ENST00000404276.6) %>% 
  summarise(counts=sum(Counts)/total.counts) %>% arrange(desc(counts))
total.counts <- sum(hadley.dat$Total.counts)
h.annotated <- hadley.dat %>% filter(Events !="")%>% group_by(Events,r.ENST00000404276.6) %>% 
  summarise(counts=sum(Total.counts)/total.counts)%>% arrange(desc(counts))

merged <- merge(j.annotated, h.annotated, by='r.ENST00000404276.6', all =T, suffixes=c(".jessie", ".hadley"))
merged <- merged %>% mutate(first_base=stringr::str_extract(r.ENST00000404276.6,"^[0-9]+"), 
                  first_base=ifelse(r.ENST00000404276.6 =="",1,first_base),
                  events = ifelse(is.na(Events.hadley), Events.jessie, Events.hadley)) %>% 
  arrange(as.numeric(first_base)) %>% select(r.ENST00000404276.6, events, counts.jessie, counts.hadley)


wb <- createWorkbook()
addWorksheet(wb, 'data')

hs <- createStyle(
  fgFill = "#DCE6F1", textDecoration = "bold",
  border = "Bottom"
)
writeData(wb, 'data', merged,headerStyle=hs)

saveWorkbook(wb, "git/gitHub/RNA-isoforms-chek2/data/summaryCHEK2.xlsx", overwrite = TRUE)
