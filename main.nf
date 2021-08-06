params.gene_expressions_fp = "/Users/timbarry/research_offsite/gasperini-2019/at-scale/processed/gene_expressions"
params.gRNA_counts_fp = "/Users/timbarry/research_offsite/gasperini-2019/at-scale/processed/gRNA_counts"
params.pairs_fp="/Users/timbarry/research_offsite/glmeiv/public/data_analysis/pairs_sample.rds"
params.gene_precomp_pod_size = 3
params.gRNA_precomp_pod_size = 3
params.pair_pod_size = 4


/*********************
* gene precomputations
*********************/
// Obtain the gene IDs.
process obtain_gene_id {
  output:
  stdout gene_id_ch_raw

  """
  Rscript -e 'pairs <- readRDS("$params.pairs_fp");
  gene_names <- unique(as.character(pairs[["gene_id"]]));
  cat(paste0(gene_names, collapse = "\n"))'
  """
}

// Transform gene_id_ch_raw into usable form.
gene_id_ch_precomp = gene_id_ch_raw.splitText().map{it.trim()}.collate(params.gene_precomp_pod_size).map{it.join(' ')}

// Run the gene precomputations.
process run_gene_precomp {
  echo true

  output:
  file '*.rds' into gene_precomp_ch_raw

  input:
  val gene_id from gene_id_ch_precomp

  """
  Rscript $projectDir/bin/run_precomp.R $params.gene_expressions_fp $gene_id
  """
}

// Create a map of (gene-id, file-path) pairs.
gene_precomp_ch = gene_precomp_ch_raw.flatten().map{file -> tuple(file.baseName, file)}


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
gRNA_id_ch_precomp = gRNA_id_ch_raw.splitText().map{it.trim()}.collate(params.gRNA_precomp_pod_size).map{it.join(' ')}

// Run the gRNA precomputations.
process run_gRNA_precomp {
  echo true

  output:
  file '*.rds' into gRNA_precomp_ch_raw

  input:
  val gRNA_id from gRNA_id_ch_precomp

  """
  Rscript $projectDir/bin/run_precomp.R $params.gRNA_counts_fp $gRNA_id
  """
}

// Create a map of (gRNA-id, file-path) pairs.
gRNA_precomp_ch = gRNA_precomp_ch_raw.flatten().map{file -> tuple(file.baseName, file)}


/************************
* gene-gRNA pair analyses
************************/
process obtain_pair_id {
  output:
  stdout all_par_ch_raw

  """
  Rscript -e 'pairs <- readRDS("$params.pairs_fp");
  gene_names <- as.character(pairs[["gene_id"]]);
  gRNA_names <- as.character(pairs[["gRNA_id"]]);
  cat(paste(gene_names, gRNA_names, collapse = "\n"))'
  """
}
// Look up the precomputation files for all pairs.
all_pairs_ch = all_par_ch_raw.splitText().map{it.trim()}.map{it.split(" ")}
all_pair_genes_labelled_ch = gene_precomp_ch.cross(all_pairs_ch).map{[it[1][1], it[1][0], it[0][1]]}
all_pairs_labelled_ch = gRNA_precomp_ch.cross(all_pair_genes_labelled_ch).map{[it[1][1], it[1][0], it[1][2], it[0][1]]}

all_pairs_labelled_ch.view()

// Buffer the all_pairs_labelled array
// all_pairs_labelled_buffered = all_pairs_labelled_ch.collate(params.pair_pod_size)

/*
// Run the gene-gRNA analysis
process run_gene_gRNA_analysis {
  echo true

  input:
  tuple val(gene_id), val(gRNA_id), val(gene_precomp_file), val(gRNA_precomp_file) from all_pairs_labelled_ch

  """
  echo $gene_id $gene_precomp_file $gRNA_id $gRNA_precomp_file
  """
}
*/

/****************
* collect results
*****************/
