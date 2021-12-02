// All arguments:
// 1. gene expressions and metadata: gene_odm, gene_metadata
// 2. gRNA counts and metadata: gRNA_odm, gRNA_metadata
// 3. gene-gRNA pairs: pairs
// 4. covariate matrix and offsets: covariate_matrix, m_offsets, g_offsets
// 5. family strings: m_fam_str, g_fam_str
// 6. results dir: result_dir
// 7. pod sizes
params.m_theta = "NA"
params.g_theta = "NA"
params.result_file_name = "result_glmeiv.rds"

/*********************
* gene precomputations
*********************/
// Obtain the gene IDs.
process obtain_gene_id {
  time "60s"

  output:
  stdout gene_id_ch_raw

  input:
  path pairs_fp from params.pairs

  """
  Rscript -e 'pairs <- readRDS("$pairs_fp");
  gene_names <- unique(as.character(pairs[["gene_id"]]));
  cat(paste0(gene_names, collapse = "\n"))'
  """
}

// Transform gene_id_ch_raw into usable form.
gene_id_ch_precomp = gene_id_ch_raw.splitText().map{it.trim()}.collate(params.gene_pod_size).map{it.join(' ')}

// Run the gene precomputations.
process run_gene_precomp {
  errorStrategy  { task.attempt <= 4  ? 'retry' : 'finish' }
  time { 20.s * params.gene_pod_size } // request 1 minute/gene of wall time

  output:
  file '*.rds' into gene_precomp_ch_raw

  input:
  path gene_odm_fp from params.gene_odm
  path gene_metadata_fp from params.gene_metadata
  path covariate_matrix_fp from params.covariate_matrix
  path m_offsets_fp from params.m_offsets
  val gene_id from gene_id_ch_precomp

  // args: 1. gene backing odm file
  // 2. gene metadata RDS file
  // 3. covariate matrix
  // 4. m offsets fp
  // 5. family string
  // 6. gene id strings
  """
  run_precomp.R $gene_odm_fp $gene_metadata_fp $covariate_matrix_fp $m_offsets_fp $params.m_fam_str $params.m_theta $gene_id
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

  input:
  path pairs_fp from params.pairs

  """
  Rscript -e 'pairs <- readRDS("$pairs_fp");
  gRNA_names <- unique(as.character(pairs[["gRNA_id"]]));
  cat(paste0(gRNA_names, collapse = "\n"))'
  """
}

// Transform gRNA_id_ch_raw into usable form.
gRNA_id_ch_precomp = gRNA_id_ch_raw.splitText().map{it.trim()}.collate(params.gRNA_pod_size).map{it.join(' ')}

// Run the gRNA precomputations.
process run_gRNA_precomp {
  errorStrategy  { task.attempt <= 4  ? 'retry' : 'finish' }
  time { 20.s * params.gRNA_pod_size } // request 1 minute/gRNA of wall time

  output:
  file '*.rds' into gRNA_precomp_ch_raw

  input:
  path gRNA_odm_fp from params.gRNA_odm
  path gRNA_metadata_fp from params.gRNA_metadata
  path covariate_matrix_fp from params.covariate_matrix
  path g_offsets_fp from params.g_offsets
  val gRNA_id from gRNA_id_ch_precomp

  // args: 1. gRNA backing odm file
  // 2. gRNA metadata RDS file
  // 3. covariate matrix
  // 4. g offsets fp
  // 5. family string
  // 6. gRNA id strings
  """
  run_precomp.R $gRNA_odm_fp $gRNA_metadata_fp $covariate_matrix_fp $g_offsets_fp $params.g_fam_str $params.g_theta $gRNA_id
  """
}

// Create a map of (gRNA-id, file-path) pairs.
gRNA_precomp_ch = gRNA_precomp_ch_raw.flatten().map{file -> tuple(file.baseName, file)}


/*******************************
* gene-gRNA pair analysis setup
*******************************/
// Obtain the pair IDs
process obtain_pair_id {
  time "60s"

  output:
  stdout all_par_ch_raw

  input:
  path pairs_fp from params.pairs

  """
  Rscript -e 'pairs <- readRDS("$pairs_fp");
  pairs <- dplyr::arrange(pairs, gene_id);
  gene_names <- as.character(pairs[["gene_id"]]);
  gRNA_names <- as.character(pairs[["gRNA_id"]]);
  cat(paste(gene_names, gRNA_names, collapse = "\n"))'
  """
}
// Look up the precomputation files for all pairs; collate.
all_pairs_ch = all_par_ch_raw.splitText().map{it.trim()}.map{it.split(" ")}
all_pair_genes_labelled_ch = gene_precomp_ch.cross(all_pairs_ch).map{[it[1][1], it[1][0], it[0][1]]}
all_pairs_labelled_ch = gRNA_precomp_ch.cross(all_pair_genes_labelled_ch).map{[it[1][1], it[1][0], it[1][2], it[0][1]]}.collate(params.pair_pod_size)
// Order the elements of the above emission, keeping only the precomputation files
def my_spread(elem_list, j) {
  out = []
  for (elem in elem_list) {
    out.add(elem[j])
  }
  return out
}
def my_spread_str(elem_list, j) {
  l_size = elem_list.size();
  out = ""
  for (i = 0; i < l_size; i ++) {
    elem = elem_list[i]
    out += (elem[j] + (i == l_size - 1 ? "" : " "))
  }
  return out
}
all_pairs_labelled_ordered = all_pairs_labelled_ch.map{[my_spread_str(it, 0), my_spread_str(it, 1), my_spread(it, 2), my_spread(it, 3)]}


/************************
* gene-gRNA pair analysis
*************************/
process run_gene_gRNA_analysis {
  echo true
  errorStrategy  { task.attempt <= 4  ? 'retry' : 'ignore' }
  time { 2.m * params.pair_pod_size }

  output:
  file 'raw_result.rds' into raw_results_ch

  input:
  path covariate_matrix_fp from params.covariate_matrix
  path gene_odm_fp from params.gene_odm
  path gene_metadata_fp from params.gene_metadata
  path m_offsets_fp from params.m_offsets
  path gRNA_odm_fp from params.gRNA_odm
  path gRNA_metadata_fp from params.gRNA_metadata
  path g_offsets_fp from params.g_offsets
  tuple val(gene_id), val(gRNA_id), file('gene_fp'), file('gRNA_fp') from all_pairs_labelled_ordered

  //args: 1. covariate_matrix
  // 2. gene_odm
  // 3. gene_metadata
  // 4. m_offsets
  // 5. gRNA_odm
  // 6. gRNA_metadata
  // 7. g_offsets
  // 8. input (gene_precomp_1, ..., gene_precomp_n), (gRNA_precomp_1, ..., gRNA_precomp_n)
  """
  run_analysis.R $covariate_matrix_fp $gene_odm_fp $gene_metadata_fp $m_offsets_fp $gRNA_odm_fp $gRNA_metadata_fp $g_offsets_fp $gene_id $gRNA_id gene_fp* gRNA_fp*
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
  file "$params.result_file_name" into collected_results_ch

  input:
  path pairs_fp from params.pairs
  file 'raw_result' from raw_results_ch.collect()

  """
  collect_results.R $params.result_file_name $pairs_fp raw_result*
  """
}
