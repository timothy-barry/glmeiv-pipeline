library(magrittr)
library(ondisc)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
gene_expressions_fp <- args[1]
gRNA_counts_fp <- args[2]
other_args <- args[seq(3L, n_args)]

# load the gene expression and gRNA count data
gene_expressions <- read_odm(odm_fp = gene_expressions_fp)
gRNA_counts <- read_odm(odm_fp = gRNA_counts_fp)

# obtain the gene IDs, gRNA IDs, gene precomp file paths, and gRNA precomp file paths.
gene_ids <- other_args[seq(from = 1L, by = 4, to = length(other_args))]
gRNA_ids <- other_args[seq(from = 2L, by = 4, to = length(other_args))]
gene_precomp_fps <- other_args[seq(from = 3L, by = 4, to = length(other_args))]
gRNA_precomp_fps <- other_args[seq(from = 4L, by = 4, to = length(other_args))]
n_pairs <- length(gene_ids)

# Loop over the pairs, loading the data and running the computation on each.
out <- lapply(seq(1, n_pairs), function(i) {
  gene_id <- gene_ids[i]; gRNA_id <- gRNA_ids[i]
  gene_precomp <- readRDS(gene_precomp_fps[i]); gRNA_precomp <- readRDS(gRNA_precomp_fps[i])
  m <- as.numeric(gene_expressions[[gene_id,]]); g <- as.numeric(gRNA_counts[[gRNA_id,]])
  data.frame(gene_id = gene_id, gRNA_id = gRNA_id, sum_m = sum(m), sum_g = sum(g))
}) %>% do.call(args = ., what = rbind)

saveRDS(out, "raw_result.rds")
