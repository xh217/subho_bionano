## Annotate complex event blocks
complex_block_annotation <- function(block,window_size=10000,USCS_gene="TxDb.Hsapiens.UCSC.hg38.knownGene") {
  library(tidyr)
  library(GenomicRanges)
  library("GenomicFeatures")
  library(USCS_gene)
  library(org.Hs.eg.db)
  txdb <- USCS_gene
  txdb.gene <-genes(txdb, columns=c("TXCHROM","TXSTART", "TXEND","GENEID","TXNAME"))
  blocks <- as.data.frame(block)
  colnames(blocks) <- "blocks"
  blocks_tmp <- separate(blocks, blocks, c("chr", "start"))
  blocks_tmp$chr <- gsub("C", "chr", blocks_tmp$chr)
  blocks_tmp$start <- as.numeric(blocks_tmp$start)
  blocks_tmp$start <- blocks_tmp$start * window_size
  blocks_tmp$end <- blocks_tmp$start + window_size
  blocks_tmp_GR <- makeGRangesFromDataFrame(blocks_tmp,
                                            keep.extra.columns = TRUE,
                                            ignore.strand = FALSE)
  hits <- findOverlaps(txdb.gene, blocks_tmp_GR)
  txdb.gene_hits <- txdb.gene[queryHits(hits), ]
  ENTREZID_gene <- AnnotationDbi::select(
                   org.Hs.eg.db,
                   keys = as.character(txdb.gene_hits$GENEID),
                   columns = c("SYMBOL")
                   )
  block_gene <- unique(ENTREZID_gene$SYMBOL)
  return (block_gene)
}
