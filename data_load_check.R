processed_dir <- paste0(.get_config_path("LOCAL_GASPERINI_2019_DATA_DIR"), "at-scale/processed/")
library(ondisc)

gene_odm_fp <- paste0(processed_dir, "gasp_scale_gene_expressions.odm")
gene_metadata_fp <- paste0(processed_dir, "gasp_scale_gene_metadata.rds")
gene_expression_odm <- read_odm(odm_fp = gene_odm_fp,
                                metadata_fp = gene_metadata_fp)

gRNA_odm_fp <- paste0(processed_dir, "gasp_scale_gRNA_counts.odm")
gRNA_metadata_fp <- paste0(processed_dir, "gasp_scale_gRNA_metadata.rds")
gRNA_counts_odm <- read_odm(odm_fp = gRNA_odm_fp,
                            metadata_fp = gRNA_metadata_fp)
