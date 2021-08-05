args <- commandArgs(trailingOnly = TRUE)
gene_expressions_fp <- args[1]
to_save_fp <- args[2]
gene_id <- args[3]

# print(args[3])

# run precomputation (in this case, simply saving the vector of expressions)
odm <- ondisc::read_odm(gene_expressions_fp)
r <- as.numeric(odm[[gene_id,]])
saveRDS(object = r, file = to_save_fp)