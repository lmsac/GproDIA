library(ggrepel)

calculate.intensity.rank = function(report_matrix, intensity_column_suffix = '.mzML',
                                    min_frequency = 0.6) {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  run_names = names(report)[endsWith(names(report), intensity_column_suffix)]
  
  missing = apply(report[, run_names], 1, function(x) sum(is.na(x))) / length(run_names)
  report = subset(report, missing <= 1 - min_frequency)
  
  result = report[, !colnames(report) %in% run_names]
  result$intensity = rowMeans(report[, run_names], na.rm = TRUE)
  result$rank = nrow(result) + 1 - rank(result$intensity, ties.method = 'first')
  
  result
}


plot.intensity.rank = function(report_matrix, ..., 
                               return_data = FALSE,
                               color_column = 'key', color_values = c('#339dff', '#ff3355'),
                               size_column = 'key', size_values = c(0.5, 2),
                               shape_column = 'key', shape_values = c(1, 18),
                               label_column = 'label') {
  intensity_rank = calculate.intensity.rank(
    report_matrix,
    ...
  )
  
  pl = ggplot(intensity_rank) + 
    geom_point(
      aes(
        x = rank, y = log10(intensity), 
        color = if (color_column %in% colnames(intensity_rank)) intensity_rank[[color_column]] else 'None',
        size = if (size_column %in% colnames(intensity_rank)) intensity_rank[[size_column]] else 'None',
        shape = if (shape_column %in% colnames(intensity_rank)) intensity_rank[[shape_column]] else 'None'
      ),
      alpha = 0.75
    ) 
  
  if (label_column %in% colnames(intensity_rank)) {
    pl = pl +
      geom_text_repel(
        data = subset(intensity_rank, !is.na(intensity_rank[[label_column]])),
        aes_string(
          x = 'rank', y = 'log10(intensity)',
          label = label_column
        ),
        size = 2, nudge_y = sample(c(-0.5, 0.5), sum(!is.na(intensity_rank[[label_column]])), replace = TRUE),
        segment.size = 0.2, min.segment.length = 0, 
        segment.color = 'darkgray', segment.alpha = 0.75
      )
  }
  
  pl = pl +
    scale_x_continuous(
      name = 'Rank'
    ) +
    scale_y_continuous(
      name = 'Log10 Intensity'
    ) +
    scale_color_manual(values = color_values) +
    scale_size_manual(values = size_values) +
    scale_shape_manual(values = shape_values) +
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
    list(plot = pl, data = intensity_rank)
  }
  else {
    pl
  }
}




report_matrices$precursor = local({
  report = report_matrices$precursor
  PRM_result = PRM_result_site_pro
  
  idx = match(
    with(
      report_matrices$precursor, 
      paste(FullPeptideName, GlycanComposition, Charge)
    ),
    with(
      PRM_result,
      paste(FullPeptideName, GlycanComposition, Charge)
    )
  )
  
  report$PRM_supported = ifelse(
    is.na(idx), 'None', 
    ifelse(PRM_result$Observed[idx], 'PRM', 'no PRM')
  )
  
  report$DIA_unique = ifelse(
    is.na(idx), 'None', 
    ifelse(
      PRM_result$DIA_only_pro[idx],
      'protein',
      ifelse(
        PRM_result$DIA_only_site[idx],
        'site',
        'siteglycan'
      )
    )
  )
  
  report
})

plot.intensity.rank(
  report_matrices$precursor[order(report_matrices$precursor$DIA_unique, decreasing = TRUE), ],
  color_column = 'DIA_unique', color_values = c('None' = '#339dff', 'siteglycan' = '#65c3ba', 'site' = '#ffa447', 'protein' = '#ff3355'),
  size_column = 'PRM_supported', size_values = c('None' = 0.05, 'no PRM' = 0.05, 'PRM' = 2),
  shape_column = 'PRM_supported', shape_values = c('None' = 20, 'no PRM' = 20, 'PRM' = 18)
)

ggsave('rank_intensity_precursor.svg', width = 6, height = 6, unit = 'cm')


report_matrices$siteglycan = local({
  report = report_matrices$siteglycan
  PRM_result = PRM_result_site_pro[order(PRM_result_site_pro$Observed, decreasing = TRUE), ]
  PRM_result = subset(PRM_result, !duplicated(paste(Proteins, ProSite, GlycanComposition)))
  
  idx = match(
    with(
      report_matrices$siteglycan, 
      paste(ProteinName, ProteinGlycoSite, GlycanComposition)
    ),
    with(
      PRM_result,
      paste(Proteins, ProSite, GlycanComposition)
    )
  )
  
  report$PRM_supported = ifelse(
    is.na(idx), 'None', 
    ifelse(PRM_result$Observed[idx], 'PRM', 'no PRM')
  )
  
  report$DIA_unique = ifelse(
    is.na(idx), 'None', 
    ifelse(
      PRM_result$DIA_only_pro[idx],
      'protein',
      ifelse(
        PRM_result$DIA_only_site[idx],
        'site',
        'siteglycan'
      )
    )
  )
  
  report
})
  
plot.intensity.rank(
  report_matrices$siteglycan[order(report_matrices$siteglycan$DIA_unique, decreasing = TRUE), ],
  color_column = 'DIA_unique', color_values = c('None' = '#339dff', 'siteglycan' = '#65c3ba', 'site' = '#ffa447', 'protein' = '#ff3355'),
  size_column = 'PRM_supported', size_values = c('None' = 0.05, 'no PRM' = 0.05, 'PRM' = 2),
  shape_column = 'PRM_supported', shape_values = c('None' = 20, 'no PRM' = 20, 'PRM' = 18)
)

ggsave('rank_intensity_siteglycan.svg', width = 6, height = 6, unit = 'cm') 


report_matrices$glycosite = local({
  report = report_matrices$glycosite
  PRM_result = PRM_result_site_pro[order(PRM_result_site_pro$Observed, decreasing = TRUE), ]
  PRM_result = subset(PRM_result, !duplicated(paste(Proteins, ProSite)))
  
  idx = match(
    with(
      report, 
      paste(ProteinName, ProteinGlycoSite)
    ),
    with(
      PRM_result,
      paste(Proteins, ProSite)
    )
  )
  
  report$PRM_supported = ifelse(
    is.na(idx), 'None', 
    ifelse(PRM_result$Observed[idx], 'PRM', 'no PRM')
  )
  
  report$DIA_unique = ifelse(
    is.na(idx), 'None', 
    ifelse(
      PRM_result$DIA_only_pro[idx],
      'protein',
      ifelse(
        PRM_result$DIA_only_site[idx],
        'site',
        'None'
      )
    )
  )
  
  report$label = ifelse(
    report$PRM_supported == 'PRM' & report$DIA_unique != 'None', 
    with(report, paste0(
      'N', ProteinGlycoSite, '@', 
      sapply(strsplit(ProteinName, '\\|'), function(s) s[2])
    )),
    NA
  )
    
  report
})


plot.intensity.rank(
  report_matrices$glycosite[order(report_matrices$glycosite$DIA_unique, decreasing = TRUE), ],
  color_column = 'DIA_unique', color_values = c('None' = '#339dff', 'site' = '#ffa447', 'protein' = '#ff3355'),
  size_column = 'PRM_supported', size_values = c('None' = 0.05, 'no PRM' = 0.05, 'PRM' = 2),
  shape_column = 'PRM_supported', shape_values = c('None' = 20, 'no PRM' = 20, 'PRM' = 18)
)

ggsave('rank_intensity_glycosite.svg', width = 12, height = 6, unit = 'cm')

