library(rtracklayer)
library(GenomicRanges)
library(gtrellis)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

path = "git/gitHub/RNA-isoforms-chek2/data/bedGraph/"

seq_depth<-import.bedGraph(paste0(path,"JessieRun10_chek2.bedgraph"))
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- 'UCSC'
CHEK2_region <- range(transcripts(edb, filter = ~symbol == "CHEK2")) +500
CHEK2_txs <- exonsBy(edb, filter = ~ symbol == 'CHEK2')

n=length(CHEK2_txs)
width=4
height=3
jes.track.height=1
had.track.height=1
lab.size=8
######### Jessie's data ######### ######### 
png("Jessie_CHEK2coverage_naturalScale.png", width=width, height=height, units = 'in',res=300)
gtrellis_layout(
  data = CHEK2_region,
  track_ylim = c(c(0, 2.2e5), c(0.5, n+0.5)),
  track_axis = c(TRUE< FALSE),
  n_track = 2,
  track_ylab = c("Coverage", "CHEK2"),
  track_height = unit.c(unit(1, 'null'), unit(jes.track.height, 'in')),
  add_name_track = TRUE,
  asist_ticks = F,
  lab_fontsize = lab.size
  )
add_lines_track(
  seq_depth, 
  seq_depth$score, 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)

add_track(panel_fun = function(gr) {
  tr = CHEK2_txs[1:23, ] # all transcripts for this gene
  for(i in seq_along(tr)) {
    # for each transcript
    current_tr_start = min(start(tr[[i]]))
    current_tr_end = max(end(tr[[i]]))
    grid.lines(c(current_tr_start, current_tr_end), c(n - i + 1, n - i + 1), 
               default.units = "native", gp = gpar(col = "#CCCCCC"))
    grid.rect(start(tr[[i]]), n - i + 1, width(tr[[i]]), 0.6,
              default.units = "native", just = "left", 
              gp = gpar(fill = "orange", col = "orange"))
  }
})
dev.off()




png("Jessie_CHEK2coverage_log10.png", width=width, height=height, units = 'in',res=300)
gtrellis_layout(
  data = CHEK2_region,
  track_ylim = c(c(0, log10(2.2e5)), c(0.5, n+0.5)),
  track_axis = c(TRUE< FALSE),
  n_track = 2,
  track_ylab = c("Coverage (log10)", "CHEK2"),
  track_height = unit.c(unit(1, 'null'), unit(jes.track.height, 'in')),
  add_name_track = TRUE,
  asist_ticks = F,
  lab_fontsize = lab.size
)
add_lines_track(
  seq_depth, 
  log10(seq_depth$score+1), 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)

add_track(panel_fun = function(gr) {
  tr = CHEK2_txs[1:23, ] # all transcripts for this gene
  for(i in seq_along(tr)) {
    # for each transcript
    current_tr_start = min(start(tr[[i]]))
    current_tr_end = max(end(tr[[i]]))
    grid.lines(c(current_tr_start, current_tr_end), c(n - i + 1, n - i + 1), 
               default.units = "native", gp = gpar(col = "#CCCCCC"))
    grid.rect(start(tr[[i]]), n - i + 1, width(tr[[i]]), 0.6,
              default.units = "native", just = "left", 
              gp = gpar(fill = "orange", col = "orange"))
  }
})
dev.off()

######### ######### ######### ######### 

######### Hadley's data ######### ######### 
b1.seq <- import.bedGraph(paste0(path,"barcode91_chek2.bedgraph"))
b2.seq <- import.bedGraph(paste0(path,"barcode92_chek2.bedgraph"))
b3.seq <- import.bedGraph(paste0(path,"barcode93_chek2.bedgraph"))
b4.seq <- import.bedGraph(paste0(path,"barcode94_chek2.bedgraph"))


png("Hadley_CHEK2coverage_naturalScale.png", width=width, height=height, units = 'in',res=300)
gtrellis_layout(
  data = CHEK2_region,
  track_ylim = c(rep(c(0, 2.2e5),4), c(0.5, n+0.5)),
  track_axis = c(TRUE< FALSE),
  n_track = 5,
  track_ylab = c(paste("Sample",1:4), "CHEK2"),
  #track_height = unit.c(unit(1, 'null'), unit(5, 'cm')),
  track_height = unit.c(unit(.25, 'null'),unit(.25, 'null'),unit(.25, 'null'),unit(.25, 'null'), unit(had.track.height, 'in')),
  asist_ticks = F,
  add_name_track = TRUE, 
  #axis_label_fontsize = lab.size,
  lab_fontsize = lab.size
)
add_lines_track(
  b1.seq, 
  b1.seq$score, 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)
add_lines_track(
  b2.seq, 
  b2.seq$score, 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)
add_lines_track(
  b3.seq, 
  b3.seq$score, 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)
add_lines_track(
  b4.seq, 
  b4.seq$score, 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)

add_track(panel_fun = function(gr) {
  tr = CHEK2_txs[1:23, ] # all transcripts for this gene
  for(i in seq_along(tr)) {
    # for each transcript
    current_tr_start = min(start(tr[[i]]))
    current_tr_end = max(end(tr[[i]]))
    grid.lines(c(current_tr_start, current_tr_end), c(n - i + 1, n - i + 1), 
               default.units = "native", gp = gpar(col = "#CCCCCC"))
    grid.rect(start(tr[[i]]), n - i + 1, width(tr[[i]]), 0.6,
              default.units = "native", just = "left", 
              gp = gpar(fill = "orange", col = "orange"))
  }
})
dev.off()

png("Hadley_CHEK2coverage_logScale.png", width=width, height=height, units = 'in',res=300)
gtrellis_layout(
  data = CHEK2_region,
  track_ylim = c(rep(c(0, log10(2.2e5)),4), c(0.5, n+0.5)),
  track_axis = c(TRUE< FALSE),
  n_track = 5,
  track_ylab = c(paste("Sample",1:4), "CHEK2"),
  #track_height = unit.c(unit(1, 'null'), unit(5, 'cm')),
  track_height = unit.c(unit(.25, 'null'),unit(.25, 'null'),unit(.25, 'null'),unit(.25, 'null'), unit(had.track.height, 'in')),
  asist_ticks = F,
  add_name_track = TRUE,
  lab_fontsize = lab.size
)
add_lines_track(
  b1.seq, 
  log10(b1.seq$score+1), 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)
add_lines_track(
  b2.seq, 
  log10(b2.seq$score+1), 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)
add_lines_track(
  b3.seq, 
  log10(b3.seq$score+1), 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)
add_lines_track(
  b4.seq, 
  log10(b4.seq$score+1), 
  area = TRUE,
  gp = gpar(fill = "gray", col = NA)
)

add_track(panel_fun = function(gr) {
  tr = CHEK2_txs[1:23, ] # all transcripts for this gene
  for(i in seq_along(tr)) {
    # for each transcript
    current_tr_start = min(start(tr[[i]]))
    current_tr_end = max(end(tr[[i]]))
    grid.lines(c(current_tr_start, current_tr_end), c(n - i + 1, n - i + 1), 
               default.units = "native", gp = gpar(col = "#CCCCCC"))
    grid.rect(start(tr[[i]]), n - i + 1, width(tr[[i]]), 0.6,
              default.units = "native", just = "left", 
              gp = gpar(fill = "orange", col = "orange"))
  }
})
dev.off()

