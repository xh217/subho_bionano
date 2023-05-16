# Retrieve xmap data with blocks
overlap_bin_mapID <- function(genome_bin, xmap, window_size = 10000, ...) {
    library(GenomicRanges)
    library(data.table)
    genome_bin <-
      read.table(
        genome_bin,
        sep = "\t",
        header = F,
        stringsAsFactors = F
      )
    colnames(genome_bin) <- c("chr", "start", "end")
    if (exists("xmap", mode = "list")) {
      contig_xmap <- xmap
    } else if (file.exists(xmap)) {
      contig_xmap <-
        utils::read.table(
          xmap,
          sep = "\t",
          header = F,
          stringsAsFactors = F
        )
    }
    colnames(contig_xmap) <- c("mapID", "chr", "start_q", "end_q", "start", "end", "strand")
    contig_xmap <- contig_xmap[, c("chr", "start", "end", "mapID")]
    contig_xmap$chr <- paste0("chr", contig_xmap$chr)
    contig_xmap[which(contig_xmap$chr == "chr23"), "chr"] = "chrX"
    ##make overlap
    refGR <-
      makeGRangesFromDataFrame(genome_bin,
                               keep.extra.columns = TRUE,
                               ignore.strand = FALSE)
    contig_xmapGR <-
      makeGRangesFromDataFrame(contig_xmap,
                               keep.extra.columns = TRUE,
                               ignore.strand = FALSE)
    hits <- findOverlaps(refGR, contig_xmapGR)
    genome_bin_tmp <- contig_xmap[subjectHits(hits),]
    contig_xmap_tmp <- genome_bin[queryHits(hits),]
    contig_xmap_tmp$chr <- gsub("chr", "C", contig_xmap_tmp$chr)
    contig_xmap_tmp$combine <-
      paste(contig_xmap_tmp$chr,
            ceiling(contig_xmap_tmp$start / window_size),
            sep = "_")
    
    combine_tmp <- as.data.frame(cbind(contig_xmap_tmp$combine, genome_bin_tmp$mapID))
    combine_tmp <- unique(combine_tmp)
    setDT(combine_tmp)[, `:=`(type_merge, paste0(as.character(V2), collapse = ",")), by = .(V1)]
    combine_tmp <- unique(as.data.frame(combine_tmp)[, c("V1", "type_merge")])
    combine_tmp$V1 <- as.character(combine_tmp$V1)
    return(combine_tmp)
  }
