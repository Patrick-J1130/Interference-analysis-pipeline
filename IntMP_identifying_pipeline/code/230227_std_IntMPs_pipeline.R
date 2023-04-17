# Library ############
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source("../code/230227_fun_library.R")
  print_time = format(Sys.time(),'%y%m%d')
  dir.create(paste0('../data/',print_time))
}

# Input ############
{
  LC_type = "LC_method_1"
  sample_info = read.xlsx("../dependence/sample_info_230227.xlsx")
  sample_info_name = sample_info %>% separate_rows(sample_name,sep = "; ") %>% distinct(sample_name) %>% pull()
  known_library = read.xlsx(paste0("../dependence/fragment_library_",LC_type,"_230227.xlsx"))
  pos_library = read_csv(paste0("../dependence/pos_transition_queue_Hilicon_6500_",LC_type,"_230227.csv"),
                         col_names = c('accession')) %>%
    left_join(known_library %>% filter(CE > 0) %>% dplyr::select(accession,rt,Q1,Q3,CE) %>% distinct())
  neg_library = read_csv(paste0("../dependence/neg_transition_queue_Hilicon_6500_",LC_type,"_230227.csv"),
                  col_names = c('accession')) %>%
    left_join(known_library %>% filter(CE < 0) %>% dplyr::select(accession,rt,Q1,Q3,CE) %>% distinct())
  # loading mzML database & peaks finding #
  mzML_files <- dir(path = paste0('../raw_data/220911_',LC_type,'_std_R/'), 
                    full.names = TRUE, recursive = F, pattern = ".mzML")
}
# Obtain chrom interf table and detected expected peak table ############
chrom_table = extract_chrom_table(mzML_files)
all_peak_list = pick_all_peaks(chrom_table,sample_info_name)
peak_table_merged = merge_peak_table(all_peak_list,sample_info,sample_type = "std",name_of_peak_finding_funs = c('centwave','matchfilter'))
peak_ls_with_chrom = add_chrom_to_peak_ls(peak_table_merged,all_peak_list)
chrom_interf = calculate_trans_sim_trans_ratio(peak_table_merged,peak_ls_with_chrom)
chrom_interf = calculate_interf_ratio(chrom_interf,peak_table_merged,peak_ls_with_chrom)
# Output chrom interf table and detected expected peak table ############
{
  write.csv(chrom_interf,paste0('../data/',print_time,'/',LC_type,'_std_chrom_interf.csv'))
  write.csv(peak_table_merged %>% filter(expected_peak == T),paste0('../data/',print_time,'/',LC_type,'_std_expected_peak_table.csv'))
}
# IntMPs identification ############
## Input ############
{
  expected_peaks = read.csv(paste0('../data/',print_time,'/',LC_type,'_std_expected_peak_table.csv'))
  chrom_interf_raw = read_csv(paste0('../data/',print_time,'/',LC_type,'_std_chrom_interf.csv'))
  pred_interf_raw = read_csv(paste0("../data/interf_prediction_",print_time,".csv"))
}
## Self data cleaning ####
{
  two_uM_std = c("HMDB0000157","HMDB0000121","HMDB0000054","HMDB0000468","HMDB0000159")
  chrom_interf = chrom_interf_raw %>%
    mutate(transition_ratio = apply(.,1,function(x){
      if(x[["interf_accession"]] %in% two_uM_std == T){
        x[["transition_ratio"]] = as.numeric(x[["transition_ratio"]])/5
      }else{
        x[["transition_ratio"]] = as.numeric(x[["transition_ratio"]])
      }
    })) 
}
## Interfering result identification ############
{
  interf_results = get_interf_results(chrom_interf = chrom_interf,pred_interf = pred_interf_raw,
                                             transition_similarity_threshold = 0.7,
                                             transition_ratio_threshold = 0.001,
                                             chrom_explained_threshold = 0.8)
}
## General_interf_table - all possible interference regardless co-elution with std 
{
  general_interf_table = get_general_interf_table(interf_results)
  general_interf_table_distinct = general_interf_table %>%
    distinct(interf_peak_id, interf_type, .keep_all=T)
  general_interf_table_summary = tabyl(general_interf_table_distinct, 
                                       interf_type)
  general_interf_table_output = general_interf_table %>%
    dplyr::select(
      interf_type,anchor_Q1,anchor_Q3,
      anchor_accession, anchor_name, 
      interf_accession, interf_name, 
      source_area, source_rt, interf_area, interf_rt, 
      transition_similarity,transition_ratio) 
  write_csv(general_interf_table_output, 
                paste0('../data/',print_time,'/',LC_type,'_general_interf_table.csv'))
    }
## Specific_interf_table - LC specific interference 
{
  specific_interf_table = get_specific_interf_table(general_interf_table,negligible_threshold = 0.005,severe_threshold = 0.5,anchor_interf_rt_diff_threshold = 0.5)
  specific_interf_table_distinct = specific_interf_table %>%
    distinct(interf_peak_id, interf_type, .keep_all=T)
  specific_interf_table_summary = tabyl(specific_interf_table_distinct, 
                                        interf_level, interf_type) 
  specific_interf_table_output = specific_interf_table %>%
    dplyr::select(
      interf_type,anchor_Q1,anchor_Q3,
      anchor_accession, anchor_name, 
      interf_accession, interf_name, 
      source_area, source_rt, interf_area, interf_rt, 
      transition_similarity,transition_ratio,normalized_overlap, interf_ratio,interf_level,anchor_rt
    )
  write_csv(specific_interf_table_output, 
            paste0('../data/',print_time,'/',LC_type,'_specific_interf_table.csv'))
  }
## Save analysis file ####
{
    saveRDS(list(general_interf_table,
                 specific_interf_table,
                 expected_peaks = expected_peaks), 
            paste0('../data/',print_time,'/',LC_type,'_interf_analysis.rds'))
  }
  

























































































