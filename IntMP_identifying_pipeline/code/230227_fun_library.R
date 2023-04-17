# Library ####
{
  library(xcms)
  library(mzR)
  library(tidyverse)
  library(openxlsx)
  library(zoo)
  library(janitor)
  print_time = format(Sys.Date(), "%y%m%d")
}
# MS2 preparation function ############
{
  # This function takes a SMILES string as input and attempts to convert it to a chemical formula.
  # It returns the molecular formula if the conversion is successful, and returns NA otherwise.
  my_SMILES2formula = function(SMILES){
    SDF = try(ChemmineR::smiles2sdf(SMILES), silent=T)
    if(inherits(SDF, "try-error")){
      warning(paste0("Invalid SMILES ", SMILES))
      return(NA)
    }
    # Try to convert the SMILES string to an SDF object
    SDF_formula = try(unname(MF(SDF,addH=T)), silent=T)
    if(inherits(SDF_formula, "try-error")){
      warning(paste0("Invalid SMILES ", SMILES))
      return(NA)
    }
    SDF_formula
  }
  
}

# predicted_interference function ####
{
  #This function calculates the number of atoms of each queried element in a 
  #given chemical formula.
  elem_num_query2 = function(formula, elem_query) 
  {
    if (!is.character(formula) | is.na(formula)) {
      return(NA)
    }
    formula <- gsub("D", "[2]H", formula)
    ende2 <- nchar(formula)
    element2 <- c()
    number2 <- c()
    j <- c(1)
    while (j <= ende2) {
      if (substr(formula, j, j) == c("[")) {
        b <- j
        while (any(substr(formula, j, j) == c("]")) != TRUE) {
          j <- c(j + 1)
        }
        k <- j
        while (any(substr(formula, j, j) == c("-", ".", 
                                              "0", "1", "2", "3", "4", "5", "6", "7", "8", 
                                              "9")) != TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        element2 <- c(element2, substr(formula, b, m))
      }
      if (any(substr(formula, j, j) == c("-", ".", "0", "1", 
                                         "2", "3", "4", "5", "6", "7", "8", "9")) != TRUE) {
        k <- j
        while (any(substr(formula, j, j) == c("-", ".", 
                                              "0", "1", "2", "3", "4", "5", "6", "7", "8", 
                                              "9")) != TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        element2 <- c(element2, substr(formula, k, m))
      }
      if (any(substr(formula, j, j) == c("-", ".", "0", "1", 
                                         "2", "3", "4", "5", "6", "7", "8", "9")) == TRUE) {
        k <- j
        while (any(substr(formula, j, j) == c("-", ".", 
                                              "0", "1", "2", "3", "4", "5", "6", "7", "8", 
                                              "9")) == TRUE) {
          j <- c(j + 1)
        }
        m <- c(j - 1)
        j <- c(j - 1)
        number2 <- c(number2, as.numeric(substr(formula, 
                                                k, m)))
      }
      j <- j + 1
    }
    
    num_query = numeric(length(elem_query))
    for(i in 1:length(element2)){
      num_query[which(elem_query==element2[i])]=number2[i]
    }
    
    return(num_query)
  }
  #This function checks if any strings in a vector are present in a 
  #comma-separated string
  fun2 = function(x,y){
    x =  unlist(strsplit(x,split = ","))
    result = any(y %in% x)
    return(result)
  }
  #This function calculates a logical matrix indicating whether each standard 
  #contains at least one atom larger of each queried element than other standard 
  calc_element_table = function(list_of_std, 
                                elem = c('C','H','O','N','P','S','Cl','I','Se','F')){
    
    elem_table_ls = list()
    for(i in 1:length(elem)){
      temp_elem_vec = list_of_std[[paste0(elem[i], '_number')]]
      elem_table_ls[[i]] = outer(temp_elem_vec,temp_elem_vec,`-`) >= 0
    }
    Reduce("&", elem_table_ls)
  }
  
  merge_tables = function(col_name = 'Q1_type', 
                          table_names = c()){
    if(length(table_names) == 0){
      return(NULL)
    }
    temp_ls = list()
    for(i in 1:length(table_names)){
      temp_table = eval(as.symbol(table_names[i]))
      temp = which(temp_table, arr.ind = T) %>%
        as.data.frame() %>%
        mutate(V3 = table_names[i])
      
      temp_ls[[i]] = temp
    }
    result = bind_rows(temp_ls) 
    colnames(result)[3] = col_name
    
    return(result)
  }
  #Add all information according to local metabolites library
  add_library_info = function(Q1Q3_table, database = new_data_base){
    database = database %>%
      dplyr::select(accession, name, formula, Q1, Q3) 
    database_anchor = database %>% rename_all(function(x){
      paste0('anchor_', x)
    })
    database_interfering = database %>% rename_all(function(x){
      paste0('interf_', x)
    })
    
    Q1Q3_table = Q1Q3_table %>%
      ungroup() %>%
      mutate(anchor_accession = list_of_std$accession[col],
             interf_accession = list_of_std$accession[row],
             .keep = 'unused') %>%
      left_join(database_anchor) %>%
      left_join(database_interfering) %>%
      select(interf_type,Q1_type,Q3_type,anchor_Q1,anchor_Q3,
             anchor_accession, anchor_name, anchor_formula,
             interf_accession, interf_name, interf_formula)
  }
}

# Peak table function ####
{
  #fuzzy_group assigns each value in a vector 
  #to a group based on its proximity to other values####
  fuzzy_group = function(data, abs_tol = 0.002, rel_tol = 10e-6){
    if(length(data) == 1){
      return(0)
    }
    # Get the indices that would sort the vector
    data_order = order(data)
    # Sort the vector
    data_sort = data[data_order]
    count = 1
    data_group = rep(1,(length(data_sort))) # Initialize group assignments for each value
    for(i in 2:length(data_sort)){
      temp_tol = max(abs_tol, data_sort[i] * rel_tol)
      if(data_sort[i]-data_sort[i-1]>temp_tol){
        count = count+1 # If the distance between adjacent values exceeds the tolerance, start a new group
      }
      data_group[data_order[i]]=count # Assign the current value to the current group
    }
    data_group
  }
  
  # generate_scantable_list reads in an mzML file and extracts the intensity values and retention times for each scan.####
  generate_scan_table = function(mzML){
      fns = sub(".mzML","", basename(mzML)) #Extract the sample name from the file name
      pd <- data.frame(sample_name = sub(basename(mzML), pattern = ".mzML",
                                         replacement = "", fixed = TRUE),
                       stringsAsFactors = FALSE) #Create a data frame to hold sample metadata
      raw_data <- readMSData(files = mzML, pdata = new("NAnnotatedDataFrame", pd),
                             mode = "onDisk") #Read in the raw data in on-disk mode
      
      spec_all = xcms::spectra(raw_data)
      #Extract scan meta data from the raw data
      scanData = raw_data@featureData@data
      pos_scan = scanData %>%
        filter(polarity == 1) %>%
        filter(grepl('experiment\\=[1|2]', spectrumId)) %>%
        pull(acquisitionNum)
      
      neg_scan = scanData %>%
        filter(polarity == 0) %>%
        filter(grepl('experiment\\=[1|2]', spectrumId)) %>%
        pull(acquisitionNum)
      
      inten_list = sapply(spec_all, function(x){
        x@intensity
      })
      rt_list = scanData$retentionTime
      
      pos = bind_rows(inten_list[pos_scan])  %>%
        rename_with(function(x){rt_list[pos_scan]}) %>%
        bind_cols(pos_library) %>%
        gather(key = 'rt', value = 'inten', -Q1, -Q3, -rt, -CE, -accession) %>%
        mutate(rt = as.numeric(rt))
      
      neg = bind_rows(inten_list[neg_scan])  %>%
        rename_with(function(x){rt_list[neg_scan]}) %>%
        bind_cols(neg_library) %>%
        gather(key = 'rt', value = 'inten', -Q1, -Q3, -rt, -CE, -accession) %>%
        mutate(rt = as.numeric(rt))
      #Merge the positive and negative polarity data frames
      merge_table = bind_rows(pos, neg)%>%
        mutate(rt= round(rt/60, 3)) %>%
        mutate(sample_name = fns)
    }
  }
  
  # Find peaks with CentWave ####
  generate_peaksWithCentWave = function(scan_table_ls){
    centwave_ls = list()
    
    for(i in c(0.4,0.6,0.8)){
      for(j in seq(0.2,0.5,0.1)){
        temp_centwave_ls = lapply(scan_table_ls, function(x){
          temp = xcms::peaksWithCentWave(x$inten, x$rt, peakwidth = c(i, i+j), 
                                         snthresh = 1, prefilter = c(2,2000),
                                         fitgauss = F,
                                         # extendLengthMSW = T,
                                         noise = 0, verboseColumns = T) %>%
            as.data.frame() %>%
            mutate(i = i,
                   j = i+j)
        })
        centwave_ls[[length(centwave_ls)+1]]=temp_centwave_ls
      }
    }
    centwave_ls
    
  }
  
  # Find peaks with MatchedFilter ####
  generate_peaksWithMatchedFilter = function(scan_table_ls){
    matchfilter_ls = list()
    for(i in seq(0.4,0.9,0.1)){
      temp_matchfilter_ls = lapply(scan_table_ls, function(x){
        temp = xcms::peaksWithMatchedFilter(x$inten, x$rt, 
                                            fwhm = i,
                                            max = 7,
                                            snthresh = 1) %>% 
          as.data.frame() %>%
          mutate(fwhm = i)
        
      })
      matchfilter_ls[[length(matchfilter_ls)+1]] = temp_matchfilter_ls
    }
    matchfilter_ls
  }
  
  # my_peak_picking funciton####
  smooth_scan_table = function(scan_table)
  {
    scan_gap = scan_table$rt[2] - scan_table$rt[1] # Calculate the time difference between two consecutive scans
    inten_noise_threshold = 2000 # Set a threshold for the intensity value below which any data point will be considered as noise and its intensity value will be set to zero
    chrom_peak_in_min = 0.5 # Specify the width of the peak in minutes
    smoothing_scan_num = floor(chrom_peak_in_min /scan_gap) # Calculate the number of scans to be used for smoothing the data
    allow_missing_in_consecutive_scans = c(missing=2,consecutive=5)
    baseline_scans = 3
    half_integration_scans = 7 ## Specify the number of scans to be used for integration. The total number of scans used for integration will be 2*half_integration_scans+1
    
    scan_table = scan_table %>%
      mutate(new_intensity = ifelse(`inten`<= inten_noise_threshold, 0,`inten`),
             rollmean = rollmean(new_intensity, smoothing_scan_num,
                                 fill = 0, align = 'center')) # Apply a rolling mean function to the new_intensity column with the specified number of scans
    return(scan_table)
    
  }
  
  
  # generate_chrom_fun ####
  generate_chrom_fun = function(peak_table, spline_n = 30){
    colnames(peak_table) = c('rt','inten')
    chrom_peak = spline(peak_table$rt,peak_table$inten,
                        n = spline_n*(length(peak_table$rt)), 
                        method = 'natural') 
    chrom_peak[['y']][chrom_peak[['y']]<0] = 0
    
    area_fun = approxfun(chrom_peak[['x']],chrom_peak[['y']])
  }
  #calculate cos overlaps between two peaks
  calculate_cos_overlaps = function(peak1, peak2, rt_expand = 0.5){
    rt_min = min(peak1$rt, peak2$rt)
    rt_max = max(peak1$rt, peak2$rt)
    
    rt_seq = seq(rt_min-rt_expand, rt_max+rt_expand, 0.001)
    
    peak1_chrom = peak1$chrom_fun(rt_seq)
    peak2_chrom = peak2$chrom_fun(rt_seq)
    peak1_chrom[is.na(peak1_chrom)] = 0
    peak2_chrom[is.na(peak2_chrom)] = 0
    
    cos_overlap_rate = peak1_chrom%*%peak2_chrom / (sqrt(peak1_chrom%*%peak1_chrom) * sqrt(peak2_chrom%*%peak2_chrom))
    if(is.na(cos_overlap_rate) == TRUE){
      cos_overlap_rate = 0
    }
    return(cos_overlap_rate)
  }
  
  # get_peaks_from_peak_group_id ####
  get_peaks_from_peak_group_id = function(peak_group_id, peak_ls){
    peak_id = sapply(peak_ls, function(x){
      if(is.na(x$peak_group_id)){
        FALSE
      } else {
        x$peak_group_id == peak_group_id
      }
    })
    peak_ls[peak_id]
  }

# main code new functions ####
{
  extract_chrom_table = function(mzML_files){
    all_peaks_table = list()
    for(mzML in mzML_files){
      peak_table = generate_scan_table(mzML)
      fns = sub(".mzML","", basename(mzML))
      all_peaks_table[[fns]] = list(
        peak_table = peak_table
      )
      
      
    }
    chrom_table = lapply(all_peaks_table, function(x){
      x$peak_table
    }) %>% bind_rows(.id='sample_name')
  }
  pick_all_peaks = function(chrom_table,sample_info_name){
    all_peak_tables = list()
    print(Sys.time())
    for(x in sample_info_name){
      tictoc::tic()
      if(x %in% chrom_table$sample_name){
        scan_table_ls = chrom_table %>%
          filter(sample_name == x) %>%
          dplyr::select(-sample_name) %>%
          split(paste(.$accession,.$Q1,.$Q3))
        # Self_developed
        {
          my_scantable = lapply(scan_table_ls, smooth_scan_table) %>% 
            bind_rows()
          }
        # peaksWithCentWave
        {
          centwave_ls = generate_peaksWithCentWave(scan_table_ls)
          
          centwave = lapply(centwave_ls, bind_rows, .id='id') %>%
            bind_rows() %>%
            separate(id, c('accession','Q1','Q3'), sep=" ") %>%
            type.convert() %>%
            mutate(origin = 'peaksWithCentWave') %>%
            group_by(accession, Q1, Q3) %>%
            mutate(rt_group = fuzzy_group(rt, 10/60)) %>%
            distinct(accession, Q1, Q3, rt_group, .keep_all=T) %>%
            dplyr::select(-rt_group) %>%
            ungroup() %>%
            arrange(accession, Q1, Q3, i, j) %>%
            filter(T)
        }
        
        # peaksWithMatchedFilter ####
        {
          matchfilter_ls = generate_peaksWithMatchedFilter(scan_table_ls)
          
          matchfilter = lapply(matchfilter_ls, bind_rows, .id='id') %>%
            bind_rows() %>%
            separate(id, c('accession','Q1','Q3'), sep=" ") %>%
            type.convert() %>%
            mutate(origin = 'peaksWithMatchedFilter') %>%
            group_by(accession, Q1, Q3) %>%
            mutate(rt_group = fuzzy_group(rt, 10/60)) %>%
            distinct(accession, Q1, Q3, rt_group, .keep_all=T) %>%
            dplyr::select(-rt_group) %>%
            ungroup() %>%
            arrange(accession, Q1, Q3, fwhm) %>%
            filter(T)
        }
        all_peak_tables[[x]] = list(
          my_scantable = my_scantable,
          centwave = centwave,
          matchfilter = matchfilter
        )
        print(x)
        tictoc::toc()
      }else{
        stop('Please check whether the names in sample_info and sample names are the same!')
      }
    }
    return(all_peak_tables)
  }
  merge_peak_table = function(all_peak_list,sample_info,sample_type = "std",name_of_peak_finding_funs = c('centwave','matchfilter')){
    peak_table_merge = bind_rows(
      lapply(all_peak_list, function(x){
        x$centwave 
      }) %>% bind_rows(.id='sample_name') %>%
        dplyr::select(accession, sample_name, rt, rtmin, rtmax, into) %>%
        mutate(algorithm = 'centwave'),
      lapply(all_peak_list, function(x){
        x$matchfilter 
      }) %>% bind_rows(.id='sample_name') %>%
        dplyr::select(accession, sample_name, rt, rtmin, rtmax, into) %>%
        mutate(algorithm = 'matchfilter')
    ) %>%
      group_by(accession, sample_name) %>%
      mutate(rt_group = fuzzy_group(rt, abs_tol = 10/60)) %>%
      group_by(accession, sample_name, rt_group) %>%
      mutate(rt = median(rt)) %>%
      ungroup() %>%
      arrange(-into) %>%
      distinct(accession, algorithm, sample_name, rt_group, .keep_all=T) %>%
      dplyr::select(accession, sample_name, into, algorithm, rt) %>%
      spread(key = 'algorithm', value = 'into') %>%
      mutate(peak_id = 1:nrow(.)) %>%
      group_by(peak_id) %>%
      mutate(inten = max(centwave,matchfilter,na.rm = T)) %>%
      ungroup()
    expected_peak_table_merge = peak_table_merge %>%
      full_join(known_library %>% 
                  filter(!is.na(rt)) %>%
                  dplyr::rename(std_rt = rt) %>%
                  dplyr::select(accession,std_rt)
      ) %>%
      filter(abs(rt - std_rt) < 0.5) %>%
      mutate(temp_sum = rowSums(.[,name_of_peak_finding_funs], na.rm=T)) %>%
      arrange(-temp_sum) %>%
      distinct(accession, sample_name, .keep_all=T) %>%
      dplyr::select(-temp_sum) %>%
      arrange(accession) %>%
      mutate(expected_peak = T) %>%
      filter(T)
    
    if(sample_type %in% c("std","standard","Std","Standard") == T){
      expected_peak_table_merge = expected_peak_table_merge %>%
        full_join(sample_info %>% filter(!is.na(sample_name)) %>% 
                    dplyr::rename(std_category = sample_name)) %>%
        separate_rows(std_category, sep='; ') %>%
        filter(sample_name == std_category)
    }
    
    peak_table_merge = left_join(peak_table_merge,expected_peak_table_merge %>% dplyr::select(peak_id,expected_peak))
    return(peak_table_merge)
  }
  add_chrom_to_peak_ls = function(peak_table_merge,all_peak_list){
    peak_ls = apply(peak_table_merge, 1, function(x){
      list(
        peak_id = as.numeric(x['peak_id']),
        peak_group_id = NA,
        accession = unname(x['accession']),
        sample_name = unname(x['sample_name']),
        rt = as.numeric(x['rt']),
        centwave = as.numeric(x['centwave']),
        matchfilter = as.numeric(x['matchfilter']),
        inten = max(as.numeric(x['centwave']), as.numeric(x['matchfilter']), na.rm=T)
      )
    })
    raw_ls = lapply(all_peak_list, function(x){
      x$my_scantable 
    }) %>% 
      bind_rows(.id='sample_name') %>%
      split(~ accession) %>%
      lapply(split, ~ sample_name)
    
    for(i in 1:length(peak_ls)){
      temp_peak = peak_ls[[i]]
      
      peak_table = raw_ls[[temp_peak$accession]][[temp_peak$sample_name]] %>%
        filter(abs(rt - temp_peak$rt) < 0.5)
      
      peak_table = peak_table[,c('rt','new_intensity')]
      
      temp_peak[['chrom_fun']] = generate_chrom_fun(peak_table, spline_n = 30)
      
      rt_seq = seq(temp_peak$rt-0.5, temp_peak$rt+0.5, 0.001)
      temp_peak[['area']] = temp_peak[['chrom_fun']](rt_seq) %>%
        sum(na.rm=T) * 0.001
      
      peak_ls[[i]] = temp_peak
    }
    return(peak_ls)
  }
  
  calculate_trans_sim_trans_ratio = function(peak_table_merge,peak_ls){
    expected_peak_table_merge = peak_table_merge %>% filter(expected_peak == T)
    cos_sim_ls = list()
    expected_peaks_ls = peak_ls[expected_peak_table_merge$peak_id]
    for(i in 1:length(expected_peaks_ls)){
      # print(i)
      temp_peak = expected_peaks_ls[[i]]
      temp_ids = sapply(peak_ls, function(x){
        x$sample_name == temp_peak$sample_name &
          abs(x$rt - temp_peak$rt) < 2
      })
      candidate_peak_ls = peak_ls[temp_ids]
      transition_similarity = sapply(candidate_peak_ls, calculate_cos_overlaps, temp_peak)
      candidate_peak_ls = candidate_peak_ls[transition_similarity > 0.2]
      if(length(candidate_peak_ls)==0){
        next
      }
      candidate_info = lapply(candidate_peak_ls, '[', c('peak_id','accession','rt', 'area')) %>%
        bind_rows() %>%
        dplyr::rename(
          interf_peak_id = peak_id,
          interf_rt = rt,
          interf_area = area,
          anchor_accession = accession
        ) %>%
        mutate(transition_similarity = transition_similarity[transition_similarity > 0.2])
      temp_chrom_interf = data.frame(
        source_peak_id = temp_peak$peak_id,
        interf_accession = temp_peak$accession,
        sample_name = temp_peak$sample_name,
        source_rt = temp_peak$rt,
        source_area = temp_peak$area
      ) %>%
        bind_cols(candidate_info) %>%
        mutate(transition_ratio = interf_area/temp_peak$area)
      cos_sim_ls[[i]] = temp_chrom_interf
    }
    
    chrom_interf = bind_rows(cos_sim_ls) %>%
      dplyr::select(sample_name, anchor_accession, interf_peak_id, interf_area, interf_rt,
                    interf_accession, source_peak_id, source_area, source_rt, 
                    transition_similarity, transition_ratio)
    return(chrom_interf)
  }
  calculate_interf_ratio = function(chrom_interf_pre,peak_table_merge,peak_ls){
    expected_peak_table_merge = peak_table_merge %>% filter(expected_peak == T)
    std_peak_ls = peak_ls[expected_peak_table_merge$peak_id]
    chrom_interf_ls = chrom_interf_pre %>%
      dplyr::select(interf_peak_id, anchor_accession) %>%
      distinct() %>%
      split(1:nrow(.))
    for(i in 1:length(chrom_interf_ls)){
      temp_id = chrom_interf_ls[[i]]$interf_peak_id
      temp_accession = chrom_interf_ls[[i]]$anchor_accession
      temp_peak = peak_ls[[temp_id]]
      temp_std_peaks_id = sapply(std_peak_ls, function(x){
        x$accession == temp_accession
      })
      temp_std_peaks = std_peak_ls[temp_std_peaks_id]
      if(length(temp_std_peaks)==0){
        temp_overlap = 1
        interf_ratio = Inf
      } else {
        temp_overlap = sapply(temp_std_peaks, calculate_cos_overlaps, temp_peak)
        temp_overlap_id = which.max(temp_overlap)
        temp_overlap = temp_overlap[temp_overlap_id]
        interf_ratio = temp_peak$area / temp_std_peaks[[temp_overlap_id]]$area
      }
      chrom_interf_ls[[i]]['normalized_overlap']=temp_overlap
      chrom_interf_ls[[i]]['interf_ratio']=interf_ratio
    }
    
    interf_table = bind_rows(chrom_interf_ls)
    chrom_interf = chrom_interf_pre %>%
      left_join(interf_table)
    return(chrom_interf)
  }
  get_interf_results = function(chrom_interf,pred_interf,transition_similarity_threshold = 0.7,transition_ratio_threshold = 0.001,chrom_explained_threshold = 0.8){
    chrom_interf = chrom_interf %>%
      filter(transition_similarity > transition_similarity_threshold) %>%
      filter(transition_ratio > transition_ratio_threshold) %>%
      left_join(sample_info %>%
                  separate_rows(sample_name,sep='; ') %>%
                  dplyr::select(accession, sample_name) %>%
                  dplyr::rename(anchor_sample_category = sample_name), 
                by=c('anchor_accession'='accession')) %>%
      left_join(sample_info %>%
                  dplyr::select(accession, sample_name) %>%
                  dplyr::rename(interf_sample_category = sample_name), 
                by=c('interf_accession'='accession')) %>%
      filter(T)
    pred_interf = pred_interf %>%
      group_by(anchor_accession, interf_accession) %>%
      mutate(predicted_pair_id = cur_group_id()) 
    
    interf_results = full_join(chrom_interf, pred_interf) %>%
      filter(!is.na(interf_peak_id)) %>% 
      mutate(pred_explained = !is.na(predicted_pair_id),
             chrom_explained = transition_similarity > chrom_explained_threshold) %>%
      mutate(chrom_pred_score = transition_similarity + 0.1*pred_explained) %>%
      arrange(-chrom_pred_score) %>%
      group_by(interf_peak_id) %>%
      mutate(rank = 1:n()) %>%
      ungroup() %>%
      filter(T)
  }
  get_general_interf_table = function(interf_results){
    general_interf_table_id = interf_results %>%
      filter(pred_explained, chrom_explained) %>%
      distinct(interf_peak_id, .keep_all=T) %>%
      filter(!mapply(grepl,anchor_sample_category,interf_sample_category)) %>%
      dplyr::select('interf_peak_id','predicted_pair_id') %>%
      filter(T)
    
    general_interf_table = interf_results %>%
      right_join(general_interf_table_id) %>%
      distinct(interf_peak_id,interf_type, .keep_all=T)
    
    return(general_interf_table)
  }  
  get_specific_interf_table = function(general_interf_table,negligible_threshold = 0.005,severe_threshold = 0.5,anchor_interf_rt_diff_threshold = 0.5){
    specific_interf_table = general_interf_table %>%
      left_join(expected_peaks %>% 
                  dplyr::select(accession,sample_name,rt) %>% 
                  dplyr::rename(anchor_accession = accession,anchor_sample_category = sample_name,anchor_rt = rt),
                by = c("anchor_accession","anchor_sample_category")) %>%
      mutate(interf_level = case_when(
        interf_ratio < negligible_threshold ~ 'negligible',
        interf_ratio >= negligible_threshold & interf_ratio < severe_threshold ~ 'slight',
        interf_ratio >= severe_threshold ~ 'severe'
      )) %>%
      #提取anchor rt
      filter(abs(anchor_rt - interf_rt) < anchor_interf_rt_diff_threshold & interf_ratio != "Inf") %>%
      filter(interf_level != 'negligible')
    return(specific_interf_table)
  }
  identify_interf_in_sample = function(chrom_interf_sample,ref_interf_table,
                                       transition_similarity_threshold = 0.8,
                                       interf_ratio_sample_min_threshold = 0.1,
                                       interf_ratio_sample_max_threshold = 10,
                                       sample_technical_or_biological_duplication = T,
                                       sample_repetitions = 3){
    sample_interf = chrom_interf_sample %>%
      left_join(ref_interf_table_distinct, 
                by = c('anchor_accession','interf_accession'),
                suffix = c('','_ref')) %>%
      filter(transition_similarity > transition_similarity_threshold) %>%
      mutate(interf_ratio_sample = transition_ratio_ref/transition_ratio) %>%
      filter(interf_ratio_sample<interf_ratio_sample_max_threshold, interf_ratio_sample>interf_ratio_sample_min_threshold)
    if(sample_technical_or_biological_duplication == T){
      sample_interf = sample_interf %>%
        group_by(interf_accession, anchor_accession) %>%
        mutate(rt_group = fuzzy_group(interf_rt, abs_tol = 15/60)) %>%
        group_by(anchor_accession,interf_accession, rt_group) %>%
        filter(sd(source_area)/mean(source_area)<0.5) %>%
        filter(length(unique(sample_name))>=sample_repetitions) %>%
        filter(T)
      
      # conslidate into peak group
      sample_interf = sample_interf %>%
        select(-source_peak_id, -interf_peak_id) %>%
        group_by(anchor_accession, interf_accession, rt_group) %>%
        mutate_at(vars(interf_area, interf_rt, source_area, source_rt,
                       transition_similarity, transition_ratio, interf_ratio_sample),
                  median) %>%
        mutate(sample_name = paste(unique(sort(sample_name)), collapse = ', ')) %>%
        distinct() %>%
        ungroup() %>%
        select(-rt_group) %>%
        filter(T) 
    }
    return(sample_interf)
  }
}


