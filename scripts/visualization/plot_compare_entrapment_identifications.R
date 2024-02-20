library(ggplot2)


plot.compare.entrapment.identifications = function(reports, level = 'precursor', 
                                           ...,
                                           total_entrapment = TRUE,
                                           return_data = FALSE) {
  entrapment_identifications = lapply(reports, function(report) {
    get.entrapment.identifications(report, level, ...)
  })
  
  data = do.call(rbind, lapply(names(entrapment_identifications), function(name) {
    data = entrapment_identifications[[name]]
    data$run = as.factor(data$run)
    data$dataset = factor(name, levels = names(entrapment_identifications))
    data
  }))
  if (total_entrapment) {
    data$entrapment = ifelse(data$entrapment == 'Non-Entrap', 'Non-Entrap', 'Entrap')
    data = aggregate(count ~ ., data, sum)
  }
  else {
    data = subset(data, count != 0)
  }
  
  pl = ggplot(
    data, 
    aes(x = run, y = ifelse(entrapment == 'Non-Entrap', count, -count))
  ) +
    geom_bar(
      aes(fill = entrapment),
      stat = 'identity', # position = position_dodge(),
      alpha = 0.6, color = NA
    ) +
    geom_text(
      aes(
        label = count, 
        y = ifelse(entrapment == 'Non-Entrap', count, - count), 
        color = entrapment,
        hjust = ifelse(xor(entrapment == 'Non-Entrap', abs(count) < 10), 'right', 'left')
      ),
      # position = position_dodge(width = 0.9),
      angle = 90,
      size = 3
      #check_overlap = TRUE
    ) +
    
    scale_fill_manual(values = c(
      'Entrap Glycan' = '#339dff', 
      'Non-Entrap' = '#65c3ba', 
      'Entrap Peptide' = '#ffa447',
      'Entrap Both' = '#ff3355',
      'Entrap' = '#ff3355'
    )) +
    scale_color_manual(values = c(
      'Entrap Glycan' = '#194e7f', 
      'Non-Entrap' = '#3c756f', 
      'Entrap Peptide' = '#ffd247',
      'Entrap Both' = '#7f192a',
      'Entrap' = '#7f192a'
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
    facet_grid(cols = vars(dataset), scales = 'free_x', space = 'free_x') +
    theme(
      axis.line.x = element_line(), 
      axis.line.y = element_line(), 
      panel.background = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.title.y = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = element_text(color = 'black'),
      strip.background = element_blank(),
      strip.text = {
        if (length(unique(data$dataset)) == 1) element_blank() 
        else element_text(color = 'black')
      },
      legend.position = 'none'
    )
  
  if (return_data) {
    list(plot = pl, data = entrapment_identifications)
  }
  else {
    pl
  }
}



if (exists('reports') && exists('output_dir')) {
  entrapment_identifications = local({
    levels = c('precursor')
    
    entrapment_identifications = lapply(levels, function(level) {
      res = plot.compare.entrapment.identifications(
        reports, level, return_data = TRUE, 
        use_struct = TRUE, use_site = TRUE,
        entrapment_peptide = ' ',
        entrapment_glycan = '[XG]'
      )
      
      res$plot = res$plot + theme(strip.text = element_blank())
      
      ggsave(
        paste0(output_dir, '/', 'bar_compare_entrapment_identification_', level, '.svg'), 
        res$plot, 
        width = 12, height = 6, unit = 'cm'
      )
      
      res$data
    })
    names(entrapment_identifications) = levels
    entrapment_identifications
  })
}

