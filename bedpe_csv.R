## convert bedpe file to csv file format
bedpe_csv <- function(x) {
  if (!file.exists(x)) {
    cat("no input file exists...")
    stop()
  }
  else{
    cat("find input file\n")
  }
  
  #read header
  cat("header for vcf ...\n")
  con <- file(x)
  open(con)
  no.body <- FALSE
  x.header <- NULL
  x.colnames <- NULL
  
 # loop over header and assign file colnames
  while (TRUE) {
    x.header.tmp <- readLines(con, n = 1)
    
    if (0 == length(x.header.tmp)) {
      no.body <- TRUE
      
      break
      
    }
    else if (!grepl("^#", x.header.tmp)) {
      break
      
    }
    else if (grepl("^#[Cc]", x.header.tmp)) {
      x.colnames <- unlist(strsplit(x.header.tmp, split = "\t| +"))
      
      x.colnames[1] <- gsub("^#", "", x.colnames[1])
    }
    else{
      x.header <- c(x.header, x.header.tmp)
    }
  }
  close(con)
  cat("read header done\n")
  
  ##read files
  df_input <- read.table(x, header = F, sep = '\t')
  colnames(df_input) <- x.colnames
  
  df <- data.frame(
    CHROM_A = character(0),
    START_A = numeric(0),
    END_A = numeric(0),
    CHROM_B = character(0),
    START_B = numeric(0),
    END_B = numeric(0),
    ID = character(0),
    QUAL = character(0),
    TYPE = character(0)
    )
  if (nrow(df_input) > 0) {
    cat('exsit SV mutations\n')
    all.bedpe <- data.frame(
      CHROM_A = gsub("chr", "", df_input$CHROM_A),
      START_A = df_input$START_A,
      END_A = df_input$END_A,
      CHROM_B = gsub("chr", "", df_input$CHROM_B),
      START_B = df_input$START_B,
      END_B = df_input$END_B,
      ID = df_input$ID,
      QUAL = ifelse(is.na(df_input$QUAL), '.', df_input$QUAL),
      #STRAND_A = df_input$STRAND_A,
      #STRAND_B = df_input$STRAND_B,
      TYPE = df_input$TYPE
    )
  } else {
    # assign empty bedpe
    all.bedpe <- df  
  }
  return(all.bedpe)
}  
