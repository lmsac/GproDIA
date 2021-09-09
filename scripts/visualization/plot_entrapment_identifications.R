library(ggplot2)

get.entrapment.identifications = function(report, level = 'precursor', 
                                          qvalue = 0.01, 
                                          ms2_peakgroup_qvalue = 0.01, glycopeptide_qvalue = 0.01,
                                          entrapment_peptide = 'HUMAN',
                                          entrapment_glycan = '[FAG]',
                                          ...) {
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  if ('peak_group_rank' %in% colnames(report)) {
    report = subset(report, peak_group_rank == 1)
  }
  
  run_names = sort(unique(report$filename))
  do.call(rbind, lapply(1:length(run_names), function(i) {
    
    report = subset(report, filename == run_names[i])
    report = report[order(report$m_score), ]
    report = subset(report, !duplicated(report[, get.id.columns(level, ...)]))
    
    peptide_entrap = grepl(entrapment_peptide, report$ProteinName)
    glycan_entrap = grepl(entrapment_glycan, 
                          if ('GlycanStruct' %in% colnames(report)) report$GlycanStruct 
                          else report$GlycanComposition)
    
    passed = report$m_score < qvalue
    if ('ms2_m_score' %in% colnames(report)) {
      passed = passed & (report$ms2_m_score < ms2_peakgroup_qvalue)
    }
    if ('m_score_glycopeptide_global' %in% colnames(report)) {
      passed = passed & (report$m_score_glycopeptide_global < glycopeptide_qvalue)
    }
    
    count = c(
      'Non-Entrap' = sum(passed & !peptide_entrap & !glycan_entrap),
      'Entrap Glycan' = sum(passed & !peptide_entrap & glycan_entrap),
      'Entrap Peptide' = sum(passed & peptide_entrap & !glycan_entrap),
      'Entrap Both' = sum(passed & peptide_entrap & glycan_entrap)
    )
    data.frame(
      run = i,
      filename = run_names[i],
      qvalue = qvalue,
      entrapment = names(count),
      count = as.numeric(count),
      stringsAsFactors = FALSE
    )
  }))
}


plot.entrapment.identifications = function(report, level = 'precursor', 
                                           ...,
                                           return_data = FALSE) {
  entrapment_identifications = get.entrapment.identifications(report, level, ...)
  
  data = entrapment_identifications
  data = subset(data, count != 0)
  data$run = as.factor(data$run)
  
  pl = ggplot(
    data, 
    aes(x = run, y = ifelse(entrapment == 'Non-Entrap', count, -count))
  ) +
    geom_bar(
      aes(fill = entrapment),
      stat = 'identity', position = position_dodge(),
      alpha = 0.6, color = NA
    ) +
    geom_text(
      aes(label = count, y = ifelse(entrapment == 'Non-Entrap', count + 40, - count - 25), color = entrapment),
      position = position_dodge(width = 0.9),
      angle = 75, size = 3,
      #check_overlap = TRUE,
      hjust = 0.5
    ) +
    
    scale_fill_manual(values = c(
      'Entrap Glycan' = '#339dff', 
      'Non-Entrap' = '#65c3ba', 
      'Entrap Peptide' = '#ffa447',
      'Entrap Both' = '#ff3355'
    )) +
    scale_color_manual(values = c(
      'Entrap Glycan' = '#194e7f', 
      'Non-Entrap' = '#3c756f', 
      'Entrap Peptide' = '#ffd247',
      'Entrap Both' = '#7f192a'
    )) +
    
    scale_y_continuous(
      name = list(
        'precursor' = '# Precursors',
        'glycopeptide' = '# Glycopeptides',
        'siteglycan' = '# Site-specific Glycans',
        'glycosite' = '# Protein Glycosites'
      )[[level]]
    ) +
    scale_x_discrete(
      name = 'Run'
    ) +
    theme(
      axis.line.x = element_line(), 
      axis.line.y = element_line(), 
      panel.background = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.title.y = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = element_text(color = 'black'),
      legend.position = 'none'
    )
  
  if (return_data) {
    list(plot = pl, data = entrapment_identifications)
  }
  else {
    pl
  }
}



if (exists('report') && exists('output_dir')) {
  entrapment_identifications = local({
    levels = c('precursor')
    
    entrapment_identifications = lapply(levels, function(level) {
      res = plot.entrapment.identifications(
        report, level, return_data = TRUE,
        use_struct = TRUE, use_site = TRUE
      )
      
      ggsave(
        paste0(output_dir, '/', 'bar_entrapment_identification_', level, '.svg'), 
        res$plot, 
        width = 6, height = 6, unit = 'cm'
      )
      
      res$data
    })
    names(entrapment_identifications) = levels
    entrapment_identifications
  })
}

