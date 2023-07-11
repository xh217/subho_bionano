coor_dataframe<-function(coordinates,window=10000){
matches <- regmatches(coordinates, gregexpr("\\d+", coordinates))[[1]]

# Extract chromosome number and positions for first junction block
first_chromosome <- matches[1]
first_start_position <- as.numeric(matches[2])*1000
first_end_position <- as.numeric(matches[2])*1000 + window

# Extract chromosome number and positions for second junction block
second_chromosome <- matches[3]
second_start_position <- as.numeric(matches[4])*1000
second_end_position<- as.numeric(matches[4])*1000 + window

# Construct dataframe
first_junction <- data.frame(chr = paste0("chr",first_chromosome),
                    start = first_start_position,
                    end = first_end_position
                    )
second_junction <- data.frame(chr = paste0("chr",second_chromosome),
                    start = second_start_position,
                    end = second_end_position
                    )  
res<-list(first_junction, second_junction)
                     
return(res)                                    
}                    
