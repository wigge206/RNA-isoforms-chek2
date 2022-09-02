library(rtracklayer)
library(GenomicRanges)
library(gtrellis)
library(EnsDb.Hsapiens.v86)
library(tidyverse)

seq_depth<-import.bedGraph("intersected2.bedgraph")
edb <- EnsDb.Hsapiens.v86
seqlevelsStyle(edb) <- 'UCSC'
CHEK2_region <- range(transcripts(edb, filter = ~symbol == "CHEK2")) +500
CHEK2_txs <- exonsBy(edb, filter = ~ symbol == 'CHEK2')

######### Jessie's data ######### ######### 
png("Jessie_CHEK2coverage_naturalScale.png", width=7, height=5, units = 'in',res=300)
n=length(CHEK2_txs)
gtrellis_layout(
  data = CHEK2_region,
  track_ylim = c(c(0, 2.2e5), c(0.5, n+0.5)),
  track_axis = c(TRUE< FALSE),
  n_track = 2,
  track_ylab = c("Coverage", "CHEK2"),
  track_height = unit.c(unit(1, 'null'), unit(5, 'cm')),
  add_name_track = TRUE
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




png("Jessie_CHEK2coverage_log10.png", width=7, height=5, units = 'in',res=300)
n=23
gtrellis_layout(
  data = CHEK2_region,
  track_ylim = c(c(0, log10(2.2e5)), c(0.5, n+0.5)),
  track_axis = c(TRUE< FALSE),
  n_track = 2,
  track_ylab = c("Coverage (log10)", "CHEK2"),
  track_height = unit.c(unit(1, 'null'), unit(5, 'cm')),
  add_name_track = TRUE
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
