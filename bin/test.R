#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
precomp_fp <- args[1L]
precomp <- readRDS(precomp_fp)
print(precomp$fam)
