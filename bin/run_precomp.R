#!/usr/bin/env Rscript

########################################
# 1. Load packages and command-line args
########################################
if (!("devtools" %in% rownames(installed.packages()))) install.packages("devtools", repos = "https://cloud.r-project.org")
if (!("ondisc" %in% rownames(installed.packages()))) devtools::install_github(repo = "timothy-barry/ondisc", upgrade = "never")
if (!("glmeiv" %in% rownames(installed.packages()))) devtools::install_github(repo = "timothy-barry/glmeiv", upgrade = "never")

library(ondisc)
args <- commandArgs(trailingOnly = TRUE)
n_args <- length(args)
odm_fp <- args[1L]
metadata_fp <- args[2L]
covariate_matrix_fp <- args[3L]
offsets_fp <- args[4L]
fam <- as.character(args[5L])
theta <- if (args[6L] == "NA") NA else as.integer(args[6L])
ids <- args[seq(7L, n_args)]

####################################################
# 2. Load ODM, covariate mat, and offset; set family
####################################################
odm <- read_odm(odm_fp = odm_fp, metadata_fp = metadata_fp)
covariate_matrix <- readRDS(covariate_matrix_fp)
offset <- readRDS(offsets_fp)
fam_obj <- if (fam == "poisson") {
  poisson()
} else if (fam == "nb") {
  MASS::negative.binomial(theta)
} else {
  stop("Family unrecognized. Use either 'poisson' or 'nb'.")
}

####################################################
# 3. Iterate through the IDs, running precomputation
####################################################
for (id in ids) {
  counts <- as.numeric(odm[[id,]])
  precomp <- glmeiv::run_glmeiv_precomputation(y = counts,
                                               covariate_matrix = covariate_matrix,
                                               offset = offset,
                                               fam = glmeiv::augment_family_object(fam_obj))
  precomp$id <- id
  saveRDS(precomp, paste0(id, ".rds"))
}
