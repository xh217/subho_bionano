# Function to check if two ranges are within 5kb of each other
within_5kb <- function(range1, range2) {
  abs(range1 - range2) <= 5000
}

# Function to merge overlapping ranges within a group
merge_range<-function(data){
  merge_range_all<-list()
  
  # Initialize a vector to keep track of processed rows
  processed <- logical(nrow(data))
  
  # Iterate through each row
  for (i in 1:(nrow(data))) {
    
    if (!processed[i]) {
      
      # Get the elements from the current row
      keep_range<-data[i,]
      
      # Check if other rows have overlapping range
      for (j in (i + 1):nrow(data)) {
        if (j<(nrow(data)+1) && !processed[j]) {
          follow_range<-data[j,]
          
          #check if there are any overlapping
          if (   within_5kb(keep_range$start1, follow_range$start1) &&
                 within_5kb(keep_range$end1, follow_range$end1) &&
                 within_5kb(keep_range$start2, follow_range$start2) &&
                 within_5kb(keep_range$end2, follow_range$end2) &&
                 keep_range$chr1 == follow_range$chr1 &&
                 keep_range$chr2 == follow_range$chr2 &&
                 keep_range$type == follow_range$type
          ) {
            
            keep_range<- data[i,]	
            processed[j]<-TRUE       	
          } 
          
        } 
      }
      #store the keeping rows
      if(exists("keep_range")){
        merge_range_all[[length(merge_range_all)+1]] <-keep_range } 
      
    } 
  }
  
  return (merge_range_all)
}
