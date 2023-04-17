# Library ############
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source("../code/230227_fun_library.R")
  print_time = format(Sys.time(),'%y%m%d')
  dir.create(paste0('../data/',print_time))
}
# Input ############
{
  sample_type = "293T"
  LC_type = "LC_method_1"
  sample_info = read.xlsx(paste0("../dependence/sample_info_",sample_type,"_230227.xlsx"))
  sample_info_name = sample_info %>% separate_rows(sample_name,sep = "; ") %>% distinct(sample_name) %>% pull()
  known_library = read.xlsx(paste0("../dependence/fragment_library_",LC_type,"_230227.xlsx"))
  pos_library = read_csv(paste0("../dependence/pos_transition_queue_Hilicon_6500_",LC_type,"_230227.csv"),
                         col_names = c('accession')) %>%
    left_join(known_library %>% filter(CE > 0) %>% dplyr::select(accession,rt,Q1,Q3,CE) %>% distinct())
  
  neg_library = read_csv(paste0("../dependence/neg_transition_queue_Hilicon_6500_",LC_type,"_230227.csv"),
                         col_names = c('accession')) %>%
    left_join(known_library %>% filter(CE < 0) %>% dplyr::select(accession,rt,Q1,Q3,CE) %>% distinct())
  
  # loading mzML database & peaks finding #
  mzML_files <- dir(path = paste0('../raw_data/220911_',LC_type,'_',sample_type,'_','R/'), 
                    full.names = TRUE, recursive = F, pattern = ".mzML")
  ref_interf_table_distinct = read.csv(paste0('../data/',print_time,'/',LC_type,'_specific_interf_table.csv')) %>%
    arrange(-source_area) %>%
    distinct(interf_accession, anchor_accession, .keep_all = T) 
}

chrom_table = extract_chrom_table(mzML_files)
all_peak_list = pick_all_peaks(chrom_table,sample_info_name)
#You can provide peak table merged obtained by using your own peaks picking function 
peak_table_merged = merge_peak_table(all_peak_list,sample_info,sample_type = sample_type,name_of_peak_finding_funs = c('centwave','matchfilter'))
peak_ls_with_chrom = add_chrom_to_peak_ls(peak_table_merged,all_peak_list)
chrom_interf = calculate_trans_sim_trans_ratio(peak_table_merged,peak_ls_with_chrom)
chrom_interf = calculate_interf_ratio(chrom_interf,peak_table_merged,peak_ls_with_chrom)

sample_interf = identify_interf_in_sample(chrom_interf,ref_interf_table_distinct,
                                                      transition_similarity_threshold = 0.8,
                                                      interf_ratio_sample_min_threshold = 0.1,
                                                      interf_ratio_sample_max_threshold = 10,
                                                      sample_technical_or_biological_duplication = T,
                                                      sample_repetitions = 3)
sample_interf_output = sample_interf %>% distinct(anchor_accession,interf_accession,.keep_all = T)
sample_interf_summary = tabyl(sample_interf,interf_type)
write.csv(sample_interf, 
          paste0('../data/',print_time,'/',LC_type,'_',sample_type,'_sample_interf_table.csv'))
write.csv(sample_interf_summary, 
          paste0('../data/',print_time,'/',LC_type,'_',sample_type,'_sample_interf_table_summary.csv'))











