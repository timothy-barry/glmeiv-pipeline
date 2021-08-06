/*********************
* gRNA precomputations
*********************/

// Obtain the gRNA IDs.
process obtain_gRNA_id {
  output:
  stdout gRNA_id_ch_raw

  """
  Rscript -e 'pairs <- readRDS("$params.pairs_fp");
  gRNA_names <- unique(as.character(pairs[["gRNA_id"]]));
  cat(paste0(gRNA_names, collapse = "\n"))'
  """
}

// Transform gRNA_id_ch_raw into usable form.
gRNA_id_ch = gRNA_id_ch_raw.splitText().map{it.trim()}.collate(params.gRNA_precomp_pod_size).map{it.join(' ')}

// Run the gRNA precomputations.
process run_gRNA_precomp {
  echo true

  output:
  file '*.rds' into gRNA_precomp_ch_raw

  input:
  val gRNA_id from gRNA_id_ch

  """
  Rscript $projectDir/bin/run_precomp.R $params.gRNA_counts_fp $gRNA_id
  """
}

// Create a map of (gRNA-id, file-path) pairs.
gRNA_precomp_ch = gRNA_precomp_ch_raw.flatten().map{file -> tuple(file.baseName, file)}


// gRNA_names <- as.character(pairs[["gRNA_id"]])
// cat(paste(gene_names, gRNA_names, collapse = "\n"))'
