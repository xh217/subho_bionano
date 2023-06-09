# optimise code based on vcf2bedpe function from Helena Winata https://rdrr.io/github/uclahs-cds/public-R-bedr/

adjust.coordinates <- function(df, info_tag, start, end) {
  library(parallel)
  
  # convert adjustments to numeric
  df$start <- as.integer(start)
  
  df$end <- as.integer(end)
  
  if (!info_tag %in% names(df)) {
    warning(paste0(
      'info tag ',
      info_tag,
      ' not found in VCF file. Coordinates are not adjusted'
    ))
    
    return(list('start' = df$start, 'end' = df$end))
    
  }
  df[[info_tag]] <- strsplit(df[[info_tag]], ',')
  
  df$valid_ci <-  lapply(df[[info_tag]], length) == 2
  
  df[[info_tag]][!df$valid_ci] <- 0
  
  df$ci_start <- do.call(rbind, lapply(df[[info_tag]], as.numeric))[, 1]
  
  df$ci_end <- do.call(rbind, lapply(df[[info_tag]], as.numeric))[, 2]
  
  df$start <- df$start + df$ci_start
  
  df$start <- ifelse(df$start >= 0, df$start, 0)
  
  df$end <- df$end + df$ci_end
  
  df$end <- ifelse(df$end >= 0, df$end, 0)
  
  return(list('start' = df$start, 'end' = df$end))
  
}

get.bedpe.id <- function(df) {
  name <- df$ID
  
  if ('EVENT' %in% names(df)) {
    name <- df$EVENT
    
  } else if (any(grepl('^Manta', name))) {
    name <- gsub(':.$', '', df$ID)
    
  } else if (any(grepl('_[12]$', name))) {
    name <- gsub('_[12]$', '', df$ID)
    
  }
  return(name)
  
}

get.strand <- function(df) {
  strand <- rep('+', nrow(df))
  
  if ('STRAND' %in% names(df)) {
    strand <- df$STRAND
    
  } else if (any(grepl('\\]|\\[', df$ALT))) {
    strand[grepl('\\[.*\\[', df$ALT)] <- '-'
    
  } else {
    catv('STRAND information is not recorded in INFO:STRAND or ALT field. Returning default: + ')
  }
  return(strand)
  
}

# optimise code based on vcf2bedpe function from Helena Winata

vcf_to_bedpe <-
  function(x,
           filename = NULL,
           header = FALSE,
           verbose = TRUE) {
    catv('CONVERT VCF TO BEDPE\n')
    if (!is.null(attr(x, 'vcf')) &&
        attr(x, 'vcf') && all(names(x) == c('header', 'vcf'))) {
      x <- x$vcf
      
    } else {
      catv(' * This is not an vcf!\n')
      stop()
      
    }
    # use UCSC chromosome convention and subtract one position for start
    x$CHROM <- gsub('chr', '', x$CHROM)
    
    x$POS <- x$POS - 1
    
    if (!'SVTYPE' %in% names(x)) {
      stop('SVTYPE column must exist in VCF file')
      
    }
    # define columns for bedpe
    bedpe.cols <- data.frame(
      CHROM_A = character(0),
      START_A = numeric(0),
      END_A = numeric(0),
      CHROM_B = character(0),
      START_B = numeric(0),
      END_B = numeric(0),
      ID = character(0),
      QUAL = character(0),
      SVTYPE = character(0),
      FILTER = character(0),
      BNGTYPE = character(0)
    )
    
    # Convert simple breakends with no MATEID
    simple.bp <- subset(x, x$SVTYPE != 'BND' & is.na(x$MATEID))
    
    if (nrow(simple.bp) > 0) {
      catv('PROCESSING SIMPLE BREAKENDS\n')
      coordsA <-
        adjust.coordinates(simple.bp, 'CIPOS', simple.bp$POS, simple.bp$POS)
      
      coordsB <-
        adjust.coordinates(simple.bp, 'CIEND', simple.bp$END, simple.bp$END)
      
      name <- get.bedpe.id(simple.bp)
      
      if (!'CHR2' %in% names(simple.bp)) {
        # simple breakpoints are intrachromosomal; CHR2 = CHROM
        simple.bp$CHR2 <- simple.bp$CHROM
        
      }
      
      simple.bedpe <- data.frame(
        CHROM_A = simple.bp$CHROM,
        START_A = coordsA$start,
        END_A = coordsA$end,
        CHROM_B = gsub('chr', '', simple.bp$CHR2),
        START_B = coordsB$start,
        END_B = coordsB$end,
        ID = name,
        QUAL = ifelse(is.na(simple.bp$QUAL), '.', simple.bp$QUAL),
        SVTYPE = simple.bp$SVTYPE,
        FILTER = simple.bp$FILTER,
        BNGTYPE = simple.bp$BNGTYPE
      )
      
    } else {
      # assign empty bedpe
      simple.bedpe <- bedpe.cols
      
    }
    
    # Convert SVTYPE=='BND' but no paired breakends
    bnd.bp_1 <-
      subset(x, x$SVTYPE == 'BND' & !grepl('_1$|_2$', x$MATEID))
    
    if (nrow(bnd.bp_1) > 0) {
      catv('PROCESSING BND but no paired breakends\n')
      coordsA <-
        adjust.coordinates(bnd.bp_1, 'CIPOS', bnd.bp_1$POS, bnd.bp_1$POS)
      
      coordsB <-
        adjust.coordinates(bnd.bp_1, 'CIEND', bnd.bp_1$END, bnd.bp_1$END)
      
      name <- get.bedpe.id(bnd.bp_1)
      
      if (!'CHR2' %in% names(bnd.bp_1)) {
        # simple breakpoints are intrachromosomal; CHR2 = CHROM
        bnd.bp_1$CHR2 <- bnd.bp_1$CHROM
        
      }
      
      bnd.bp_1.bedpe <- data.frame(
        CHROM_A = bnd.bp_1$CHROM,
        START_A = coordsA$start,
        END_A = coordsA$end,
        CHROM_B = gsub('chr', '', bnd.bp_1$CHR2),
        START_B = coordsB$start,
        END_B = coordsB$end,
        ID = name,
        QUAL = ifelse(is.na(bnd.bp_1$QUAL), '.', bnd.bp_1$QUAL),
        SVTYPE = bnd.bp_1$SVTYPE,
        FILTER = bnd.bp_1$FILTER,
        BNGTYPE = bnd.bp_1$BNGTYPE
      )
      
    } else {
      # assign empty bedpe
      bnd.bp_1.bedpe <- bedpe.cols
      
    }
    
    
    # Convert paired breakends with MATEID
    bnd.bp <-
      subset(x, x$SVTYPE == 'BND' & grepl('_1$|_2$', x$MATEID))
    
    if (nrow(bnd.bp) > 0) {
      catv('PROCESSING BND paired breakends with MATEID\n')
      rownames(bnd.bp) <- bnd.bp$ID
      
      if ('MATEID' %in% names(bnd.bp)) {
        bnd.bp <- subset(bnd.bp,!is.na(bnd.bp$MATEID))
        
        bnd_pair <- list()
        
        for (id in rownames(bnd.bp)) {
          if (!id %in% bnd_pair) {
            bnd_pair[[id]] <- bnd.bp[id, 'MATEID']
            
          }
        }
        var <- bnd.bp[names(bnd_pair),]
        
        mates <- bnd.bp[unlist(unname(bnd_pair)), ]
        
        if (!'STRAND' %in% names(bnd.bp) &
            'MATESTRAND' %in% names(bnd.bp)) {
          # for delly v0.7.8 because the INFO:MATESTRAND is read while the INFO:STRAND is not
          # STRAND is not in the header so it is not read by bedr
          var$STRAND <- mates$MATESTRAND
          
          mates$STRAND <- var$MATESTRAND
          
        }
        name <- get.bedpe.id(var)
        
        coordsA <-
          adjust.coordinates(var, 'CIPOS', var$POS,  var$POS)
        
        coordsB <-
          adjust.coordinates(mates, 'CIPOS', mates$POS, mates$POS)
        
        svtype <-
          ifelse(!is.null(var$SIMPLE_TYPE), var$SIMPLE_TYPE, 'BND')
        
        bnd.bedpe <- data.frame(
          CHROM_A = var$CHROM,
          START_A = coordsA$start,
          END_A = coordsA$end,
          CHROM_B = mates$CHROM,
          START_B = coordsB$start,
          END_B = coordsB$end,
          ID = name,
          QUAL = var$QUAL,
          SVTYPE = svtype,
          FILTER = var$FILTER,
          BNGTYPE = var$BNGTYPE
        )
        
      } else {
        stop('MATEID is not present in the VCF file INFO field')
        
      }
    } else {
      bnd.bedpe <- bedpe.cols
      
    }
    bedpe_df <- rbind(simple.bedpe, bnd.bp_1.bedpe, bnd.bedpe)
    
    bedpe_df <- bedpe_df[with(bedpe_df, order(CHROM_A, START_A)),]
    
    if (!is.null(filename)) {
      write.table(
        bedpe_df,
        file = filename,
        col.names = header,
        row.names = FALSE,
        sep = '\t',
        quote = FALSE
      )
      
    }
    return(bedpe_df)
    
  }

catv <- function(x) {
  verbose <- get0("verbose", parent.frame(1))
  if (!exists("verbose")) {
    verbose <- TRUE
  }
  if (verbose) {
    cat(x)
    
  }
}
