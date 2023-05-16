##Extract contig junctions

extract_junctions <- function(one_complex_event,window_size=10000,...) {
  library(GenomicRanges)
  library(dplyr)
  library(stringr)
  library(gtools)
  library(tidyr)
  colnames(one_complex_event) <-c("mapID",
                                   "chr",
                                   "start_q",
                                   "end_q",
                                   "start_r",
                                   "end_r",
                                   "strand")
  one_complex_event <- one_complex_event[order(one_complex_event$mapID,
                                               one_complex_event$chr,
                                               one_complex_event$start_q), ]
  one_complex_event$end_r <- as.numeric(one_complex_event$end_r)
  tmp <- one_complex_event
  window = window_size
  nchr <- unique(tmp$chr)
  note_ID <-c(
      "C1",
      "C2",
      "C3",
      "C4",
      "C5",
      "C6",
      "C7",
      "C8",
      "C9",
      "C10",
      "C11",
      "C12",
      "C13",
      "C14",
      "C15",
      "C16",
      "C17",
      "C18",
      "C19",
      "C20",
      "C21",
      "C22",
      "CX"
    )
  note_ID_used <- note_ID[c(nchr)]
  nchr_df = vector(mode = "list", length = length(nchr))
  nchr_all_df <- c()
  for (i in 1:length(nchr))
  {
    nchr_tmp_df <- tmp[which(tmp$chr == nchr[i]), ]
    nchr_tmp_start <- min(nchr_tmp_df$start_r)
    nchr_tmp_end <- max(nchr_tmp_df$end_r)
    nchr_tmp_set <- seq(nchr_tmp_start, nchr_tmp_end, window)
    
    # Create empty data frame
    nchr_df[[i]] <- data.frame(matrix(NA,           
                        nrow = length(nchr_tmp_set),
                        ncol = 3))
    nchr_df[[i]]$X1 <- nchr[i]
    nchr_df[[i]]$X2 <- nchr_tmp_set + 1
    nchr_df[[i]]$X3 <- nchr_tmp_set + window
    nchr_df[[i]]$X4 <- paste0(note_ID_used[i], "_", ceiling(nchr_df[[i]]$X2 / window))
    colnames(nchr_df[[i]]) <- c("chr", "start", "end", "note")
    nchr_df[[i]]$chr <- paste0("chr", nchr_df[[i]]$chr)
  }
  
  data_tmp <- do.call(rbind, nchr_df)
  
  ## for each mapID
  mapid <- unique(tmp$mapID)
  data_all <- c()
  #one mapID have one contig
  for (i in 1:length(mapid)) {
    contig <- tmp[which(tmp$mapID == mapid[i]), c("chr", "start_r", "end_r", "strand")]
    if (nrow(contig) == 1) {
      contig_start <- min(contig$start_r)
      contig_end <- max(contig$end_r)
      #window<-ceiling(tmp_end/tmp_start) ##how to set the window size??
      window = window_size
      contig_set <- seq(contig_start, contig_end, window)
      
      ##make a dataframe with a window gap for each contig
      # Create empty data frame
      contig_tmp <- data.frame(matrix(NA,    
                               nrow = length(contig_set),
                                ncol = 4))
      contig_tmp$X1 <- unique(contig$chr)
      contig_tmp$X2 <- contig_set + 1
      contig_tmp$X3 <- contig_set + window
      contig_tmp$X4 <- as.character(contig$strand)
      colnames(contig_tmp) <- c("chr", "start", "end", "strand")
      contig_tmp$chr <- paste0("chr", contig_tmp$chr)
      contig_tmp$block_rank <- 1
    }
    
    ## one mapID have multiple contigs
    if (nrow(contig) > 1) {
      contig_tmp <- c()
      for (j in 1:nrow(contig)) {
        contig_start <- min(contig[j, ]$start_r)
        contig_end <- max(contig[j, ]$end_r)
        #window<-ceiling(tmp_end/tmp_start) ##how to set the window size
        window = window_size
        contig_set <- seq(contig_start, contig_end, window)
        
        ##make a dataframe with a window gap for each contig
        contig_tmp1 <- data.frame(matrix(NA,
                                   nrow = length(contig_set),
                                   ncol = 4))
        contig_tmp1$X1 <- unique(contig[j, ]$chr)
        contig_tmp1$X2 <- contig_set + 1
        contig_tmp1$X3 <- contig_set + window
        contig_tmp1$X4 <- as.character(contig[j, ]$strand)
        colnames(contig_tmp1) <- c("chr", "start", "end", "strand")
        contig_tmp1$chr <- paste0("chr", contig_tmp1$chr)
        
        ##add contig block order
        contig_tmp1$block_rank <- j                 
        contig_tmp <- rbind(contig_tmp, contig_tmp1)
      }
    }
    ## make overlap
    refGR <- makeGRangesFromDataFrame(data_tmp,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = FALSE)
    testGR <- makeGRangesFromDataFrame(contig_tmp,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = FALSE)
    hits <- findOverlaps(refGR, testGR)
    
    ## Calculate the percentage of overlap and apply a filter to the hits
    overlaps <- pintersect(refGR[queryHits(hits)], testGR[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(testGR[subjectHits(hits)])
    hits <- hits[percentOverlap > 0.5]
    data_final <- contig_tmp[subjectHits(hits), ]
    data_final$mapID <- mapid[i]
    data_final$note <- data_tmp[queryHits(hits), ]$note
    
    ## reverse block orders if strand is "-"
    blocks <- unique(data_final$block_rank)
    all <- c()
    for (k in blocks) {
      block_k <- data_final[which(data_final$block_rank == k), ]
      if (block_k$strand == "+") {
        block_k$note <- block_k$note
      } else {
        block_k$note <- rev(block_k$note)
      }
      all <- rbind(all, block_k)
    }
    data_all <- rbind(data_all, all)
  }
  
  ## Extract continued bins with/without gaps
  un_blocks <- unique(data_all$mapID)
  all_final <- c()
  ##set gap between bins
  gap_blocks <- 5
  for (i in 1:length(un_blocks)) {
    un_contig <- data_all[which(data_all$mapID == un_blocks[i]), c("mapID", "block_rank", "note")]
    un_contig$block_rank <- as.factor(as.character(un_contig$block_rank))
    
    ## split each block with a list
    un_contig_list <- split(un_contig, un_contig$block_rank)
    
    ## one mapID has one contig
    if (length(un_contig_list) == 1) {
      df_tmp <- un_contig_list[[1]]
      df_tmp$num <- as.numeric(gsub("A", "", df_tmp$note))
      df_tmp$diff <- ave(
                         df_tmp$num,
                         factor(df_tmp$block_rank),
                         FUN = function(x)
                         c(NA, diff(x))
                         )
      block_gap <- which(df_tmp$diff > gap_blocks)
      if (length(block_gap) > 0) {
        df_block <- df_tmp$note[block_gap]
        df_block_forward <- df_tmp$note[block_gap - 1]
        df_block_forward <- paste0("-", df_block_forward, "//")                          ##add connect symbol
        df_gap_merge <- paste0(c(rbind(df_block_forward, df_block)), collapse = "")      ##gather gap blocks
        df_all <- paste0(df_tmp$note[1], df_gap_merge, "-", df_tmp$note[nrow(df_tmp)])   ##add head and tail bins
        all <- data.frame(ID = un_blocks[i], seq = df_all)
      } else{
        h1 <- df_tmp$note[1]                                                             ##extract first row for one block
        h2 <- df_tmp$note[nrow(df_tmp)]                                                  ##extract last row for one block
        block_h <- paste(h1, h2, sep = "-") #combine first and last block
        all <- data.frame(ID = un_blocks[i], seq = block_h)
      }   
    }
    else{
    
      ##mapID has multiple contigs
      combine <- c()
      for (h in 1:length(un_contig_list)) {
        n_row <- nrow(un_contig_list[[h]])
        df_tmp <- un_contig_list[[h]]
        df_tmp$num <- as.numeric(gsub("A", "", df_tmp$note))
        df_tmp$diff <- ave(
                           df_tmp$num,
                           factor(df_tmp$block_rank),
                           FUN = function(x)
                           c(NA, diff(x))
                           )
        block_gap <- which(df_tmp$diff > gap_blocks)
        if (length(block_gap) > 0) {
          df_block <- df_tmp$note[block_gap]
          df_block_forward <- df_tmp$note[block_gap - 1]
          df_block_forward <- paste0("-", df_block_forward, "//")                         ##add connect symbol
          df_gap_merge <- paste0(c(rbind(df_block_forward, df_block)), collapse = "")     ##gather gap blocks
          df_all <- paste0(df_tmp$note[1], df_gap_merge, "-", df_tmp$note[nrow(df_tmp)])  ##add head and tail bins
          combine <- paste(c(combine, df_all), collapse = "|")                            ##combine all blocks
        } else {
          n_row <- nrow(un_contig_list[[h]])
          h1 <- un_contig_list[[h]]$note[1]                                               ##extract first row for one block
          h2 <- un_contig_list[[h]]$note[n_row]                                           ##extract last row for one block
          block_h <- paste(h1, h2, sep = "-")                                             ##combine first and last block
          combine <- paste(c(combine, block_h), collapse = "|")                           ##combine all blocks
        }
        all <- data.frame(ID = un_blocks[i], seq = combine)
      }
    }
    all_final <- rbind(all_final, all)
  }
  
  tmp_seq <- sapply(all_final$seq, function(x)
      str_extract_all(x, "([A-Z][0-9]+_[0-9]+\\|[A-Z][0-9]+_[0-9]+)"))
  tmp_seq <- sapply(tmp_seq, function(x)
       paste0(x, collapse = ","))
  all_final$uncontinued_seq <- tmp_seq
  final_res <- all_final %>%
    mutate(note = strsplit(as.character(uncontinued_seq), "\\,")) %>%
    unnest(note) %>%
    as.data.frame %>%
    subset(select = c(ID, note))
  colnames(final_res) <- c("mapID", "note")
  
  order_event <- function(x) {
    event_split <- strsplit(x, split = "\\|")
    event_unlist <- unlist(event_split)
    event_sort <- mixedsort(event_unlist)
    event_reorder <- paste0(event_sort, collapse = "|")
    return(event_reorder)
  }
  
  tmp_note <- sapply(final_res$note, order_event)
  final_res$note <- as.vector(tmp_note)
  return(final_res)
}
