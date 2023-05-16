# convert bedpe to bed file format
bionano_bedpe_bed <- function(bedpe) {
  bedpe_file <- read.vcf(bedpe, split.info = TRUE)
  df_bio <- vcf_to_bedpe(bedpe_file)
  
  if (nrow(df_bio) == 0) {
    cat("error:no input file exsit\n")
  } else{
    df_frame <- data.frame(
      chr = character(0),
      start = numeric(0),
      end = numeric(0),
      type = character(0)
    )
    if (nrow(df_bio) > 0) {
      df_bio$CHROM_A <- as.character(df_bio$CHROM_A)
      df_bio$START_A <- as.numeric(as.character(df_bio$START_A))
      df_bio$END_A <- as.numeric(as.character(df_bio$END_A))
      df_bio$TYPE <- as.character(df_bio$BNGTYPE)
      df_bio$CHROM_B <- as.character(df_bio$CHROM_B)
      df_bio$START_B <- as.numeric(as.character(df_bio$START_B))
      df_bio$END_B <- as.numeric(as.character(df_bio$END_B))
      df_bio_A <- c()
      df_bio_A <- data.frame(
        chr = paste0("chr", df_bio$CHROM_A),
        start = df_bio$START_A,
        end = df_bio$END_A,
        type = df_bio$TYPE
      )
      df_bio_B <- c()
      df_bio_B <- data.frame(
        chr = paste0("chr", df_bio$CHROM_B),
        start = df_bio$START_B,
        end = df_bio$END_B,
        type = df_bio$TYPE
      )
      df_bio_all <- rbind(df_bio_A, df_bio_B)
    }
  }
}
