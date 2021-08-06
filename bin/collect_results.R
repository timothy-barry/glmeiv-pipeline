args <- commandArgs(trailingOnly = TRUE)
out <- do.call(rbind, lapply(args, function(fp) readRDS(fp)))
saveRDS(out, "result.rds")