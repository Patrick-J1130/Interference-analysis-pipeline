# Library ############
{
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  source("../code/230227_fun_library.R")
  print_time = format(Sys.time(),'%y%m%d')
  dir.create(paste0('../data/',print_time))
}
# Input ############
{
  known_library = read.xlsx("../dependence/fragment_library_LC_method_1_230227.xlsx")
}

  # List handling ############
  {
    ## elements calculating####
    {
      elems = sort(c('C','H','O','N','S','P','Cl','Se','I','F'))
      
      elem_df = sapply(known_library$formula, elem_num_query2, elems) %>%
        t() %>%
        as.data.frame() %>%
        rename_all(function(x){paste0(elems, '_number')})
      
      list_of_std = known_library %>%
        bind_cols(elem_df) %>%
        filter(T)
    }
  }
## MS2 handling ####
{
  MS2_result = read.csv(paste0("../dependence/MS2_fragment_230227.csv"))
  MS2_df = MS2_result %>%
    dplyr::select(-source) %>%
    group_by(accession,polarity) %>%
    mutate(fragment_mz = round(fragment_mz, 0)) %>%
    mutate(MS2_mass = paste(unlist(as.numeric(fragment_mz)),collapse = ",")) %>%
    dplyr::select(-fragment_mz) %>%
    ungroup()%>%
    distinct(accession,polarity,MS2_mass,.keep_all = T)
  # 将正负离子扫描模式的二级质谱图区分开
  # 将二级质谱信息整合回list_of_std 表中准备进行matrix匹配
  list_of_std = list_of_std %>%
    left_join(MS2_df %>% filter(polarity == 'Negative'), by = c("accession")) %>%
    mutate(MS2_mass = replace_na(MS2_mass, 0)) %>%
    dplyr::rename(MS2_neg_mass = MS2_mass) %>%
    dplyr::select(-polarity) %>%
    left_join(MS2_df %>% filter(polarity == 'Positive'), by = c("accession")) %>%
    mutate(MS2_mass = replace_na(MS2_mass, 0)) %>%
    dplyr::rename(MS2_pos_mass = MS2_mass) %>%
    dplyr::select(-polarity)
  
}
## ion cleaning ####
{
  list_of_std = list_of_std %>%
    mutate(Q3_neg = ifelse(sign(CE)==-1, round(Q3, digits=0), NA),
           Q3_pos = ifelse(sign(CE)==1, round(Q3, digits=0), NA)) %>%
    mutate(mass = round(mass, digits = 0))
}

write_csv(list_of_std,paste0("../data/list_of_std_",print_time,".csv"))

# Matrix ####
{
  ## Q1 MS2 fragment prerequisites 
  {
    mass = list_of_std$mass
    MS2_prerequisites_mass = outer(mass,mass,`-`) > 0
    MS2_prerequisites_element = calc_element_table(list_of_std, 
                                                   elem = c('C','H','O','N','P','S','Cl','I','Se','F'))
  }
  
  ## Q1 interfering 
  {
    # Same formula
    formula = list_of_std$formula
    Q1_same_formula = outer(formula,formula,'==') # d_formula
    # Same mz ## Should include above formula
    mass = list_of_std$mass
    Q1_same_mz = outer(mass,mass,`-`) == 0
    # Interfering ion's fragment formula ~ anchor formula
    # predicted spectra only
    # Q1_MS2_formula_pos = outer(list_of_std$MS2_pos_fragment_formula,list_of_std$formula,Vectorize(fun2)) & 
    #   MS2_prerequisites_mass & MS2_prerequisites_element
    # Q1_MS2_formula_neg = outer(list_of_std$MS2_neg_fragment_formula,list_of_std$formula,Vectorize(fun2)) & 
    #   MS2_prerequisites_mass & MS2_prerequisites_element
    # Interfering ion's fragment mass ~ anchor mass 
    # include predicted and experimental spectra, ## Predicted ones should include above formula
    Q1_MS2_mz_pos = outer(list_of_std$MS2_pos_mass,list_of_std$mass+1,Vectorize(fun2)) & 
      MS2_prerequisites_mass 
    Q1_MS2_mz_neg = outer(list_of_std$MS2_neg_mass,list_of_std$mass-1,Vectorize(fun2)) & 
      MS2_prerequisites_mass 
    # Isotope mass
    Q1_IS_mz = outer(mass,mass,`-`) == -1
    # 1. Same formula
    Q1_isomer = c('Q1_same_formula')
    # 2. Same mass
    Q1_isobar = c('Q1_same_mz')
    # 3. Fragment producing same mass/formula
    Q1_fragment = c('Q1_MS2_mz_pos','Q1_MS2_mz_neg')
    # 4. Isotope mass
    Q1_isotope = c('Q1_IS_mz')
  }
  
  ## Q3 interfering 
  {
    # Same Q3 mz
    Q3_same_mz_neg = outer(list_of_std$Q3_neg,list_of_std$Q3_neg,`-`) == 0
    Q3_same_mz_pos = outer(list_of_std$Q3_pos,list_of_std$Q3_pos,`-`) == 0
    
    # Interfering ion's fragment mz = anchor Q3 mz
    # Predicted and experimental
    Q3_MS2_neg = outer(list_of_std$MS2_neg_mass,list_of_std$Q3_neg,Vectorize(fun2))
    Q3_MS2_pos = outer(list_of_std$MS2_pos_mass,list_of_std$Q3_pos,Vectorize(fun2))  
    
    # Interfering ion Q3 mz = anchor Q3 mz + 1
    Q3_IS_mz_neg = outer(list_of_std$Q3_neg,list_of_std$Q3_neg,`-`) == -1
    Q3_IS_mz_pos = outer(list_of_std$Q3_pos,list_of_std$Q3_pos,`-`) == -1
    
    # Interfering ion's fragment mz = anchor Q3 mz + 1
    Q3_IS_MS2_neg = outer(list_of_std$MS2_neg_mass,list_of_std$Q3_neg - 1,Vectorize(fun2))
    Q3_IS_MS2_pos = outer(list_of_std$MS2_pos_mass,list_of_std$Q3_pos - 1,Vectorize(fun2)) 

    # Q3 interfering is divided into three types
    # 1. Same Q3, or interfering ion's MS2 predicted/experimental fragment = anchor's Q3
    Q3_same = c('Q3_same_mz_neg','Q3_same_mz_pos',
                'Q3_MS2_neg','Q3_MS2_pos')
    # 3. isotope producing same Q3 
    Q3_isotope = c('Q3_IS_mz_neg','Q3_IS_mz_pos',
                   'Q3_IS_MS2_neg','Q3_IS_MS2_pos')
  }
}

# Predicted interference ####
{
  # Each column is an anchor
  # Table 1 - isomer
  {
    Q1_table = merge_tables(col_name = 'Q1_type',
                            table_names = Q1_isomer)
    
    Q3_table = merge_tables(col_name = 'Q3_type',
                            table_names = c(Q3_same)) 
    
    table_isomer = full_join(Q1_table,Q3_table) %>%
      filter(complete.cases(.)) %>%
      mutate(interf_type = 'isomer') %>%
      add_library_info(database = known_library)
  }
  
  # Table 2 - isobar
  {
    Q1_table = merge_tables(col_name = 'Q1_type',
                            table_names = Q1_isobar)
    
    Q3_table = merge_tables(col_name = 'Q3_type',
                            table_names = c(Q3_same)) 
    
    table_isobar = full_join(Q1_table,Q3_table) %>%
      filter(complete.cases(.)) %>%
      mutate(interf_type = 'isobar') %>%
      add_library_info(database = known_library) %>%
      anti_join(table_isomer, by = c('anchor_accession',
                                     'interf_accession'))  
  }
  
  # Table 3 - fragment
  {
    Q1_table = merge_tables(col_name = 'Q1_type',
                            table_names = Q1_fragment)
    
    Q3_table = merge_tables(col_name = 'Q3_type',
                            table_names = c(Q3_same)) 
    
    table_fragment = full_join(Q1_table,Q3_table) %>%
      filter(complete.cases(.)) %>%
      filter((grepl('pos',Q1_type) & grepl('pos',Q3_type)) |
               (grepl('neg',Q1_type) & grepl('neg',Q3_type))) %>%
      mutate(interf_type = 'fragment') %>%
      add_library_info(database = known_library)
  }
  
  # Table 4 - isotope
  {
    Q1_table = merge_tables(col_name = 'Q1_type',
                            table_names = Q1_isotope)
    
    Q3_table = merge_tables(col_name = 'Q3_type',
                            table_names = c(Q3_same, Q3_isotope)) 
    
    table_isotope = full_join(Q1_table,Q3_table) %>%
      filter(complete.cases(.)) %>%
      mutate(interf_type = 'isotope') %>%
      add_library_info(database = known_library)
  }
  
  table_total = bind_rows(table_isomer,
                          table_isobar,
                          table_fragment,
                          table_isotope) %>%
    distinct()
  
  table_total_distinct = table_total %>% 
    distinct(anchor_accession,interf_accession,.keep_all = T) 
  
  table_total_summary = janitor::tabyl(table_total_distinct,interf_type) %>%
    slice(c(3,2,1,4))
  
}

# Output ####
{
  write_csv(table_total %>% dplyr::select(-Q1_type,-Q3_type,-anchor_formula,-interf_formula),paste0("../data/interf_prediction_",print_time,".csv"))
  write_csv(table_total_summary,paste0("../output/interf_prediction_summary_",print_time,".csv"))
}




















