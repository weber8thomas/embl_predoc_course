##Initialize parameters
library('RCircos')
library('RColorBrewer')
data(RCircos.Histogram.Data)
data(RCircos.Heatmap.Data)
data(RCircos.Link.Data)
data("UCSC.HG38.Human.CytoBandIdeogram")

chr.exclude <- "chrY";
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 0;
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);

rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.cyto <- RCircos.Get.Plot.Ideogram();
rcircos.position <- RCircos.Get.Plot.Positions();
RCircos.List.Plot.Parameters()

rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$heatmap.color <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(256);
# RCircos.Reset.Plot.Parameters(rcircos.params);
RCircos.List.Plot.Parameters();

##Plot ideogram
out.file <- "RCircosDemoHumanGenome.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot();
dev.off()
