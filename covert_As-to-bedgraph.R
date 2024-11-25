library(dplyr)

# Set working directory
setwd("~/path/to/working_directory")

# Load output data from consecutive A's Perl code 
data <- read.csv('As_20bpSliding.csv', sep='\t')

# Assign new column names
colnames(data) <- c('Position', 'As')

# Print column names to check 
print(colnames(data))

# Change cutoff if needed
prepare_bedgraph <- function(data, metric_col, output_file, cutoff = 10) { 
  bedgraph_data <- data %>%
    # Filter values based on the cutoff
    filter(!!sym(metric_col) > cutoff) %>%
    # Prepare the correct Bedgraph format
    select(chrom = Position, start = Position, end = Position, value = !!sym(metric_col))
  
  # Set chromosome name
  bedgraph_data$chrom <- "Pf3D7_07_v3"
  
  # Adjust 'start' to be 0-based and 'end' to be one base after 'start'
  bedgraph_data$start <- bedgraph_data$start - 1
  bedgraph_data$end <- bedgraph_data$start + 1
  
  # Write to file
  write.table(bedgraph_data, file = output_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)
}

# Create BEDGRAPH file with cutoff for consecutive As
prepare_bedgraph(data, "As", "consecutive_As_cutoff.bedgraph", cutoff = 10)
