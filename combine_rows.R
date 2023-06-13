# Function to combine rows based on shared elements
combine_rows <- function(data) {
  # Initialize a list to store the combined rows
  combined_rows <- list()
  
  # Initialize a vector to keep track of processed rows
  processed <- logical(nrow(data))
  
  # Iterate through each row
  for (i in 1:(nrow(data)-1)) {
  	
    if (!processed[i]) {
      # Get the elements from the current row
      elements <- unlist(strsplit(data[i,], ";"))
      
      # Initialize a vector to store the combined elements
      combined_elements <- elements
      
      # Check if other rows have shared elements
      for (j in (i + 1):nrow(data)) {
        if (!processed[j]) {
          # Get the elements from the current row being compared
          compare_elements <- unlist(strsplit(data[j, ], ";"))
          
          # Check if there are shared elements
          if (any(compare_elements %in% combined_elements)) {
            # Combine the elements
            combined_elements <- unique(c(combined_elements, compare_elements))
            
            # Mark the row as processed
            processed[j] <- TRUE
          }
        }
      }
      
      # Store the combined row
      combined_rows[[length(combined_rows) + 1]] <- paste(combined_elements, collapse = ";")
    }
  }
  
  # Return the combined rows
  return(combined_rows)
}
