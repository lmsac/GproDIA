get.id.columns = function(level, use_struct = FALSE, use_site = FALSE, include_upper_levels = FALSE, ...) {
  id_columns = list(
    'precursor' = c(
      'FullPeptideName', 
      if (use_struct) 'GlycanStruct' else 'GlycanComposition',
      if (use_site) 'GlycanSite' else NULL,
      'Charge'
    ),
    'glycopeptide' = c(
      'Sequence', 
      if (use_struct) 'GlycanStruct' else 'GlycanComposition',
      if (use_site) 'GlycanSite' else NULL
    ),
    'siteglycan' = c(
      'ProteinName',
      'ProteinGlycoSite',
      if (use_struct) 'GlycanStruct' else 'GlycanComposition'
    ),
    'glycosite' = c(
      'ProteinName',
      'ProteinGlycoSite'
    ),
    'peptide' = c(
      'Sequence'
    ),
    'protein' = c(
      'ProteinName'
    )
  )
  
  if (include_upper_levels) {
    index = match(level, names(id_columns))
    index_pep = match('peptide', names(id_columns))
    if (index < index_pep) {
      unique(unlist(id_columns[seq(index, index_pep - 1)]))
    }
    else {
      unique(unlist(id_columns[seq(index, length(id_columns))]))
    }
  }
  else {
    id_columns[[level]]
  }
}


get.identifications = function(report, level = 'precursor', ...) {
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  if ('peak_group_rank' %in% colnames(report)) {
    report = subset(report, peak_group_rank == 1)
  }
  
  run_names = sort(unique(report$filename))
  result = lapply(run_names, function(run) {
    unique(report[report$filename == run, get.id.columns(level, ...), drop = FALSE])
  })
  
  names(result) = run_names
  result
}


get.identifications.from.matrix = function(report_matrix, level = 'precursor', intensity_column_suffix = '.mzML', ...) {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  run_names = names(report)[endsWith(names(report), intensity_column_suffix)]
  result = lapply(run_names, function(run) {
    unique(report[!is.na(report[[run]]), get.id.columns(level, ...), drop = FALSE])
  })
  
  names(result) = run_names
  result
}


get.report.matrix = function(report, intensity = 'ms2') {
  if ('peak_group_rank' %in% colnames(report)) {
    report = subset(report, peak_group_rank == 1)
  }
  
  if (intensity == 'ms2') {
    intensity_column = 'Intensity'
  }
  else if (intensity == 'ms1') {
    intensity_column = 'aggr_prec_Peak_Area'
  }
  if ('decoy' %in% colnames(report)) {
    f = formula(
      decoy + Sequence + FullPeptideName + 
        GlycanStruct + GlycanComposition + GlycanSite +
        Charge + ProteinName + ProteinGlycoSite 
      ~ filename
    )
  }
  else {
    f = formula(
      Sequence + FullPeptideName + 
        GlycanStruct + GlycanComposition + GlycanSite +
        Charge + ProteinName + ProteinGlycoSite 
      ~ filename
    )
  }
  report_matrix = reshape2::dcast(
    report, 
    f, 
    value.var = intensity_column, 
    fun.aggregate = function(x) { 
      if (length(x) == 0) -1 else sum(x)
    }
  )
  
  lapply(unique(report$filename), function(run) {
    report_matrix[report_matrix[[run]] == -1, run] <<- NA
  })
  report_matrix
}


summarise.intensities = function(report_matrix, level = 'glycopeptide', n_top = 3, intensity_column_suffix = '.mzML', ...) {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  library(dplyr)
  
  id_columns = get.id.columns(level = level, ...)
  
  report = report %>%
    mutate(Intensity_mean = rowMeans(select(., ends_with(intensity_column_suffix)), na.rm = TRUE)) %>%
    group_by_at(vars(one_of(id_columns))) %>%
    mutate(!!paste0('used_for_', level, '_quantity') := row_number(desc(Intensity_mean)) <= n_top) %>%
    mutate(Intensity_mean = NULL)
  
  report = report %>% 
    filter(!!as.symbol(paste0('used_for_', level, '_quantity'))) %>%
    select(one_of(id_columns), ends_with(intensity_column_suffix)) %>%
    summarise_at(vars(ends_with(intensity_column_suffix)), function(x) {
      if (all(is.na(x))) {
        NA
      }
      else {
        sum(x, na.rm = TRUE)
      }
    }) # %>% 
    # mutate_at(vars(ends_with(intensity_column_suffix)), funs(replace(., . <= 0, NA)))
  
  report = report_matrix %>%
    select(one_of(get.id.columns(level = level, include_upper_levels = TRUE, ...))) %>%
    distinct() %>%
    right_join(report, by = id_columns)
  
  report
}


normalize.intensities.global = function(report_matrix, by = 'median', intensity_column_suffix = '.mzML') {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  if (is.character(by)) {
    by = match.fun(by)
  }
  
  quantity_columns = names(report)[endsWith(names(report), intensity_column_suffix)]
  
  ref = sapply(quantity_columns, function(x) {
    by(report[[x]], na.rm = TRUE)
  })
  ref = ref / mean(ref)
  lapply(quantity_columns, function(x) {
    report[[x]] <<- report[[x]] / ref[x]
  })
  
  report
}

normalize.intensities.quantiles = function(report_matrix, intensity_column_suffix = '.mzML') {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  quantity_columns = names(report)[endsWith(names(report), intensity_column_suffix)]
  intensity = limma::normalizeQuantiles(as.matrix(report[quantity_columns]))
  report[, quantity_columns] = intensity
  
  report
}

transform.intensities = function(report_matrix, multiplier = 1, offset = 0, intensity_column_suffix = '.mzML') {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  quantity_columns = names(report)[endsWith(names(report), intensity_column_suffix)]
  
  if (!is.null(names(multiplier))) {
    lapply(quantity_columns, function(x) {
      report[[x]] <<- report[[x]] * multiplier[x]
    })
  }
  else if (length(multiplier) == 1) {
    lapply(quantity_columns, function(x) {
      report[[x]] <<- report[[x]] * multiplier
    })
  }
  else {
    lapply(1:length(quantity_columns), function(i) {
      report[[quantity_columns[i]]] <<- report[[quantity_columns[i]]] * multiplier[i]
    })
  }
  
  if (!is.null(names(offset))) {
    lapply(quantity_columns, function(x) {
      report[[x]] <<- report[[x]] + offset[x]
    })
  }
  else if (length(offset) == 1) {
    lapply(quantity_columns, function(x) {
      report[[x]] <<- report[[x]] + offset
    })
  }
  else {
    lapply(1:length(quantity_columns), function(i) {
      report[[quantity_columns[i]]] <<- report[[quantity_columns[i]]] + offset[i]
    })
  }
  
  report
}


get.report.matrices = function(report, intensity = 'ms2',
                               normalize_by = 'median', 
                               precursor_n_top_features = 2147483647,
                               glycopeptide_n_top_precursors = 3,
                               siteglycan_n_top_glycopeptides = 3,
                               glycosite_n_top_siteglycans = 2147483647,
                               ...) {
  report_matrix = get.report.matrix(report, intensity = intensity)
  report_matrix_normalized = normalize.intensities.global(report_matrix, by = normalize_by)
  
  precursor_matrix = summarise.intensities(
    report_matrix_normalized, level = 'precursor', 
    n_top = precursor_n_top_features, 
    ...
  )
  glycopeptide_matrix = summarise.intensities(
    precursor_matrix, level = 'glycopeptide', 
    n_top = glycopeptide_n_top_precursors,
    ...
  )
  siteglycan_matrix = summarise.intensities(
    glycopeptide_matrix, level = 'siteglycan', 
    n_top = siteglycan_n_top_glycopeptides,
    ...
  )
  glycosite_matrix = summarise.intensities(
    siteglycan_matrix, level = 'glycosite', 
    n_top = glycosite_n_top_siteglycans,
    ...
  )
  
  report_matrices = list(
    precursor = precursor_matrix,
    glycopeptide = glycopeptide_matrix,
    siteglycan = siteglycan_matrix,
    glycosite = glycosite_matrix
  )
  
  lapply(names(report_matrices), function (level) {
    report_matrix = report_matrices[[level]] 
    report_matrix = subset(
      report_matrix, 
      !duplicated(report_matrix[get.id.columns(level = level)])
    )
    
    identifications = get.identifications(report, level = level)
    id_names = apply(report_matrix[colnames(identifications[[1]])], 1, paste, collapse = ' ')
    lapply(names(identifications), function(run) {
      report_matrix[!(
        id_names %in% apply(identifications[[run]], 1, paste, collapse = ' ')
      ), run] <<- NA
    })
    
    report_matrices[[level]] <<- report_matrix
  })
  
  report_matrices
}
  

calculate.cv = function(report_matrix, run_groups, intensity_column_suffix = '.mzML') {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  run_names = names(report)[endsWith(names(report), intensity_column_suffix)]
  
  result = report[, !colnames(report) %in% run_names]
  sapply(1:length(run_groups), function(i) {
    index = grep(run_groups[[i]], run_names)
    
    cv = apply(report[, run_names[index]], 1, function(x) {
      sd(x) / mean(x)
    })
    
    result[[names(run_groups)[i]]] <<- cv
  })
  
  result
}
