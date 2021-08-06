gasp_offsite <- .get_config_path("LOCAL_GASPERINI_2019_DATA_DIR")
odm_fp <- paste0(gasp_offsite, "at-scale/processed/gene_expressions")
odm <- ondisc::read_odm(odm_fp = odm_fp)
odm[[,1]]
