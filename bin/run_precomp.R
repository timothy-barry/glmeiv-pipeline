library(ondisc)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
expressions_fp <- args[1L]
metadata_fp <- args[2L]
ids <- args[seq(3L, n_args)]

odm <- read_odm(expressions_fp, metadata_fp)

# loop over gene ids, performing precomputation and saving the result
for (id in ids) {
  r <- as.numeric(odm[[id,]])
  to_save_fp <- paste0(id, ".rds")
  saveRDS(object = r, file = to_save_fp)
}
