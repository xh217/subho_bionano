# Make circos plot for block junctions
junction_circos_plot <- function(one_complex_ta_cov,link_cut_off = 4, ...) {
  library(reshape2)
  library(pals)
  library(circlize)
  library(RColorBrewer)
  
  df <- one_complex_ta_cov
  df$note <- as.character(df$note)
  df$mapID <- as.character(df$mapID)
  df_ma <- as.data.frame(reshape2::dcast(mapID ~ note, data = df, length)[, -1] > 0)
  df_ma_convert <- t(df_ma) %*% as.matrix(df_ma)
  df_ma_convert[lower.tri(df_ma_convert, diag = T)] <- 0
  df_final = data.frame( from = rep(rownames(df_ma_convert), times = ncol(df_ma_convert)),
                         to = rep(colnames(df_ma_convert), each = nrow(df_ma_convert)),
                         value = as.vector(as.matrix(df_ma_convert)),
                         stringsAsFactors = FALSE
                       )
  df_final <- df_final[which(df_final$value > 0), ]
  
  if (nrow(df_final) > 0) {
    data_cov_all <- c()
    for (i in 1:nrow(df_final)) {
      tmp <- df[which(df$note == df_final[i, 1] | df$note == df_final[i, 2]), ]
      co_mapID <- reshape2::dcast(tmp, mapID ~ note, length)
      colnames(co_mapID) <- c("mapID", "V1", "V2")
      co_mapID <- co_mapID[which(co_mapID$V1 > 0 & co_mapID$V2 > 0), ]
      co_mapID_cov <- data.frame(from = df_final[i, 1],
                                 to = df_final[i, 2],
                                 value = sum(co_mapID$V1))
      data_cov_all <- rbind(data_cov_all, co_mapID_cov)
    }
    
data_cov_all <- unique(data_cov_all[order(data_cov_all$value), ])

if (nrow(data_cov_all[which(data_cov_all$value > 0), ])) {
  data_cov_all$from <- as.character(data_cov_all$from)
  data_cov_all$to <- as.character(data_cov_all$to)
  
  #calculate how many grids
  grid_n <- unique(c(data_cov_all$from, data_cov_all$to))
  
  ##define grid colors
  grid.col = rep("grey", length(grid_n))
  names(grid.col) <- grid_n
  
  ##define link color
  num_colors <- length(unique(data_cov_all$value))
  colors <- brewer.pal(num_colors, "Blues")
  col <- colors[findInterval(data_cov_all$value, seq(
                min(data_cov_all$value),
                max(data_cov_all$value),
                length.out = num_colors))]
  
  col[data_cov_all[[3]] < link_cut_off] = "lightgrey"
    
  circos.clear()
  chordDiagram( data_cov_all,
                grid.col = grid.col,
                col = col,
                transparency = 0,
                link.zindex = rank(data_cov_all[[3]]),
                annotationTrack = "grid",
                preAllocateTracks = 1
              )
  circos.trackPlotRegion( track.index = 1,
                          panel.fun = function(x, y) {
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             sector.name = get.cell.meta.data("sector.index")
                             circos.text(
                             mean(xlim),
                             ylim[1],
                             sector.name,
                             facing = "clockwise",
                             niceFacing = TRUE,
                             adj = c(-0.3, 0),
                             cex = 0.7
                           )
    #circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
    },
    bg.border = NA
  )
  }
 } else NULL

}
    
