##Initialize parameters
library('RCircos')
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
out.file <- "RCircosDemoHumanGenome_delly.pdf";
# pdf(file=out.file, height=8, width=8, compress=TRUE);
# RCircos.Set.Plot.Area();
# RCircos.Chromosome.Ideogram.Plot();
# dev.off()

##Plot translocation using ribbon
data(RCircos.Ribbon.Data)
TL_link <- RCircos.Ribbon.Data[1,]
TL_link$chromA <- c("chr22")
TL_link$chromStartA <- c(37934425)
TL_link$chromEndA <- c(37934425)
TL_link$chromB <- c("chr13")
TL_link$chromStartB <- c(20412250)
TL_link$chromEndB <- c(20412250)

# out.file <- "RCircosDemoHumanGenome.pdf";
pdf(file=out.file, height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area();
RCircos.Chromosome.Ideogram.Plot();

track.num <- 1;
RCircos.Link.Plot(TL_link, track.num, TRUE);
RCircos.Ribbon.Plot(ribbon.data=TL_link, track.num=9, by.chromosome=TRUE, twist=FALSE);
dev.off()