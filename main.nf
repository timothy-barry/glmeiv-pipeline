params.gene_expressions_fp = "/Users/timbarry/research_offsite/gasperini-2019/at-scale/processed/gene_expressions"
params.gRNA_counts_fp = "/Users/timbarry/research_offsite/gasperini-2019/at-scale/processed/gRNA_counts"
params.pairs_fp="/Users/timbarry/research_offsite/glmeiv/public/data_analysis/pairs_sample.rds"

// Define functions to be used throughout.
def my_collate(l) {
    out = ""
    for (gene_id in l) {
      out = out + gene_id + " "
    }
    return out
}


// Obtain the gene IDs.
process obtain_gene_id {
  output:
  stdout gene_id_ch_raw

  """
  Rscript $projectDir/bin/create_pair_channel.R $params.pairs_fp
  """
}

// Transform gene_id_ch_raw into usable form.
gene_id_ch = gene_id_ch_raw.splitText().map{it.trim()}.gene_id_ch.collate(3).map{my_collate(it)}
gene_id_ch.view()

/*
process run_gene_precomp {
  echo true

  output:
  tuple val(gene_id), file('precomp.rds') into gene_precomp_ch

  input:
  val gene_id from gene_id_ch

  """
  Rscript $projectDir/bin/run_precomp.R $params.gene_expressions_fp precomp.rds $gene_id
  """
}
*/
