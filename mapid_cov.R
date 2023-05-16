## Obtain coverage information for each mapID
mapid_cov <- function(xmap, cov) {
  library(tidyr)
  library(dplyr)
  cov <- read.table(xmap,
                    header = F,
                    sep = ' ',
                    stringsAsFactors = F)
  colnames(cov) <- c("mapID", "cov")
  complex_mapID <- read.table(cov,
                    header = T,
                    sep = '\t',
                    stringsAsFactors = F)
  complex_mapID_tmp <- complex_mapID$col1
  all_complex_mapID <- c()
  for (i in 1:length(complex_mapID_tmp)) {
    note <- strsplit(complex_mapID_tmp[i], "\\,") %>%
            as.data.frame
    note$complex_index <- i
    colnames(note) <- c("mapID", "complex_index")
    all_complex_mapID <- rbind(all_complex_mapID, note)
  }
  all_complex_mapID$mapID <- as.integer(as.character(all_complex_mapID$mapID))
  complex_mapID_cov <- left_join(all_complex_mapID, cov)
  complex_mapID_cov <- complex_mapID_cov[, c("complex_index", "mapID", "cov")]
  return (complex_mapID_cov)
}
