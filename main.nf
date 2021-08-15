// All arguments:
// 1. gene expressions and metadata: gene_odm, gene_metadata
// 2. gRNA counts and metadata: gRNA_odm, gRNA_metadata
// 3. gene-gRNA pairs: pairs
// 4. covariate matrix and offsets: covariate_matrix, m_offsets, g_offsets
// 5. family strings: m_fam_str, g_fam_str
// 6. results dir: result_dir
// 7. pod sizes
// params.gene_pod_size = 3
// params.gRNA_pod_size = 3
// params.pair_pod_size = 3

/*********************
* gene precomputations
*********************/
// Obtain the gene IDs.
process obtain_gene_id {
  time "60s"

  output:
  stdout gene_id_ch_raw

  """
  Rscript -e 'pairs <- readRDS("$params.pairs");
  gene_names <- unique(as.character(pairs[["gene_id"]]));
  cat(paste0(gene_names, collapse = "\n"))'
  """
}

// Transform gene_id_ch_raw into usable form.
gene_id_ch_precomp = gene_id_ch_raw.splitText().map{it.trim()}.collate(params.gene_pod_size).map{it.join(' ')}

// Run the gene precomputations.
process run_gene_precomp {
  time { 20.s * params.gene_pod_size } // request 1 minute/gene of wall time

  output:
  file '*.rds' into gene_precomp_ch_raw

  input:
  val gene_id from gene_id_ch_precomp

  // args: 1. gene backing odm file
  // 2. gene metadata RDS file
  // 3. covariate matrix
  // 4. m offsets fp
  // 5. family string
  // 6. gene id strings
  """
  Rscript $projectDir/bin/run_precomp.R $params.gene_odm $params.gene_metadata $params.covariate_matrix $params.m_offsets $params.m_fam_str $gene_id
  """
}

// Create a map of (gene-id, file-path) pairs.
gene_precomp_ch = gene_precomp_ch_raw.flatten().map{file -> tuple(file.baseName, file)}


/*********************
* gRNA precomputations
*********************/
// Obtain the gRNA IDs.
process obtain_gRNA_id {
  time "60s"

  output:
  stdout gRNA_id_ch_raw

  """
  Rscript -e 'pairs <- readRDS("$params.pairs");
  gRNA_names <- unique(as.character(pairs[["gRNA_id"]]));
  cat(paste0(gRNA_names, collapse = "\n"))'
  """
}

// Transform gRNA_id_ch_raw into usable form.
gRNA_id_ch_precomp = gRNA_id_ch_raw.splitText().map{it.trim()}.collate(params.gRNA_pod_size).map{it.join(' ')}

// Run the gRNA precomputations.
process run_gRNA_precomp {
  time { 10.s * params.gRNA_pod_size } // request 1 minute/gRNA of wall time

  output:
  file '*.rds' into gRNA_precomp_ch_raw

  input:
  val gRNA_id from gRNA_id_ch_precomp

  // args: 1. gRNA backing odm file
  // 2. gRNA metadata RDS file
  // 3. covariate matrix
  // 4. m offsets fp
  // 5. family string
  // 6. gRNA id strings
  """
  Rscript $projectDir/bin/run_precomp.R $params.gRNA_odm $params.gRNA_metadata $params.covariate_matrix $params.m_offsets $params.g_fam_str $gRNA_id
  """
}

// Create a map of (gRNA-id, file-path) pairs.
gRNA_precomp_ch = gRNA_precomp_ch_raw.flatten().map{file -> tuple(file.baseName, file)}


/************************
* gene-gRNA pair analyses
************************/
// Obtain the pair IDs
process obtain_pair_id {
  time "60s"

  output:
  stdout all_par_ch_raw

  """
  Rscript -e 'pairs <- readRDS("$params.pairs");
  pairs <- dplyr::arrange(pairs, gene_id);
  gene_names <- as.character(pairs[["gene_id"]]);
  gRNA_names <- as.character(pairs[["gRNA_id"]]);
  cat(paste(gene_names, gRNA_names, collapse = "\n"))'
  """
}
// Look up the precomputation files for all pairs.
all_pairs_ch = all_par_ch_raw.splitText().map{it.trim()}.map{it.split(" ")}
all_pair_genes_labelled_ch = gene_precomp_ch.cross(all_pairs_ch).map{[it[1][1], it[1][0], it[0][1]]}
all_pairs_labelled_ch = gRNA_precomp_ch.cross(all_pair_genes_labelled_ch).map{[it[1][1], it[1][0], it[1][2], it[0][1]]}

// Buffer the all_pairs_labelled array
all_pairs_labelled_buffered = all_pairs_labelled_ch.collate(params.pair_pod_size).map{it.flatten()}.map{it.join(' ')}

// Run the gene-gRNA analysis
process run_gene_gRNA_analysis {
  errorStrategy  { task.attempt <= 4  ? 'retry' : 'finish' }
  time { 2.m * params.pair_pod_size }

  output:
  file 'raw_result.rds' into raw_results_ch

  input:
  val input from all_pairs_labelled_buffered

  //args: 1. covariate_matrix
  // 2. gene_odm
  // 3. gene_metadata
  // 4. m_offsets
  // 5. gRNA_odm
  // 6. gRNA_metadata
  // 7. g_offsets
  // 8. input (gene_id, gRNA_id, gene_precomp, gRNA_precomp) tuples
  """
  Rscript $projectDir/bin/run_analysis.R $params.covariate_matrix $params.gene_odm $params.gene_metadata $params.m_offsets $params.gRNA_odm $params.gRNA_metadata $params.g_offsets $input
  """
}


/****************
* collect results
*****************/
process collect_results {
  time { 5.m * task.attempt * task.attempt }
  errorStrategy 'retry'
  maxRetries 4
  publishDir params.result_dir, mode: "copy"

  output:
  file 'result_glmeiv.rds' into collected_results_ch

  input:
  file 'raw_result' from raw_results_ch.collect()

  """
  Rscript $projectDir/bin/collect_results.R $params.pairs raw_result*
  """
}
