args <- commandArgs(trailingOnly = TRUE)
# get file path
fp <- if (is.na(args[1])) "/Users/timbarry/research_offsite/glmeiv/public/data_analysis/pairs_sample.rds" else args[1] 

# load the pairs
pairs <- readRDS(fp)

# print to stdout
gene_names <- 
cat(paste0(as.character(pairs$gene_id), collapse = "\n"))
