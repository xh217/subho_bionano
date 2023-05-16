# Convert complex event to a table 
complex_event_table <- function(complex_event) {
  separator <- ","
  # Loop through the list and add the separator to each sublist
  for (i in 1:length(complex_event$contigs)) {
    complex_event$contigs[[i]] <- paste(complex_event$contigs[[i]], collapse = separator)
  }
  
  for (i in 1:length(complex_event$blocks)) {
    complex_event$blocks[[i]] <- paste(complex_event$blocks[[i]], collapse = separator)
  }
  final_cmp <- data.frame(
      mapID = do.call(rbind, complex_event$contigs),
      blocks = do.call(rbind, complex_event$blocks)
    )
}
