library(karyoploteR)
library(GenomicRanges)
library(regioneR)

repeats_cons <- toGRanges("/Users/edolzhenko/projects/2024/Q2/PlatinumPedigree/pipelines/tandem-repeats/output/consistent_trs.bed")
repeats_incons <- toGRanges("/Users/edolzhenko/projects/2024/Q2/PlatinumPedigree/pipelines/tandem-repeats/output/inconsistent_trs.bed")


pp <- getDefaultPlotParams(plot.type=2)
kp <- plotKaryotype(plot.type=2, genome = "hg38", plot.params = pp)

kpDataBackground(kp, data.panel = 1, col="#AACBFF")
kpDataBackground(kp, data.panel = 2, col="#FFAACB")

kpPlotDensity(kp, data.panel = 1, window.size = 100000, data = repeats_cons, col="red")
kpPlotDensity(kp, data.panel = 2, window.size = 100000, data = repeats_incons, col="blue")
