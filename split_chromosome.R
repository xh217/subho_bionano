# Split chromosome into 10kb windows with contextual information
 split_chromosome <- function(chromosome, start, end, window_size = 10000, ...) {
   # Calculate the number of 10kb windows in the chromosome
   n_windows <- ceiling((end - start + 1) / window_size)

   # Create a dataframe with n_windows rows
   df_new <- data.frame(matrix(ncol = 4, nrow = n_windows))
   colnames(df_new) <- c("chromosome", "start", "end", "context")

   # Populate the dataframe with the 10kb windows
   for (i in 1:n_windows) {
     window_start <- start + (i-1) * window_size
     window_end <- min(start + i * window_size - 1, end)
     df_new[i,] <- c(chromosome, window_start, window_end, ...)
   }

   return(df_new)
 }
 
