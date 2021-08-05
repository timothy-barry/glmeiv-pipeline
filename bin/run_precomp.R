args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
gene_expressions_fp <- args[1]
gene_ids <- args[seq(2L, n_args)]

odm <- ondisc::read_odm(gene_expressions_fp)

# loop over gene ids, performing precomputation and saving the result
for (gene_id in gene_ids) {
  r <- as.numeric(odm[[gene_id,]])
  to_save_fp <- paste0(gene_id, ".rds")
  saveRDS(object = r, file = to_save_fp)
}
