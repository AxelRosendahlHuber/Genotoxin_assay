# Functions accompanying the research protocol: 
# Plot structural variants

# helper function for plot_circos: 
#generate matching dataframe in bed-style format linking structural variant breakends
makeMatchingSVbeds = function(linx_list) {
  linx_list = rbindlist(linx_list)
  bed_list = list(
    bed1 = data.frame(chr = paste0("chr",linx_list$ChrStart),
                      start = linx_list$PosStart, end = linx_list$PosStart + 1, 
                      Type = linx_list$Type),
    bed2 = data.frame(chr = paste0("chr", linx_list$ChrEnd),
                      start = linx_list$PosEnd, end = linx_list$PosEnd + 1,
                      Type = linx_list$Type))
  
  return(bed_list)
}

# plot structural variants with LINX output *linx.vis_sv_data.tsv as input. 
plot_circos = function(linx_list) { 
  
  SV_beds = makeMatchingSVbeds(linx_list)
  SV_beds$color = SV_beds$bed1$Type %>% dplyr::recode(DEL = "maroon", INV = "darkblue", BND = "black", DUP = "#4daf4a")
  
  circos.initializeWithIdeogram(plotType = NULL)
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, col = "black")
    circos.text(mean(xlim), mean(ylim), 
                labels =  gsub('chr', '', chr), 
                cex = 0.7, col = "white",
                facing = "inside", niceFacing = TRUE)
  }, track.height = 0.15, bg.border = NA)
  plot = circos.genomicLink(region1 = SV_beds$bed1, SV_beds$bed2,col = SV_beds$color, lwd = 2)
  
  # plot legend
  lgd_links = Legend(title = "Structural variants", at = unique(SV_beds$bed1$Type), legend_gp = gpar(fill = unique(SV_beds$color)))
  lgd_list_horizontal = packLegend(lgd_links, direction = "horizontal")
  draw(lgd_list_horizontal, x = unit(4, "mm"), y = unit(95, "mm"), just = c("left", "top"))
}
