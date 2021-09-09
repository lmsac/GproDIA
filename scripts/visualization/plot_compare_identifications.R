library(ggplot2)

plot.compare.run.identifications = function(reports, level, ..., return_data = FALSE) {
  run_identifications = lapply(reports, function(report) get.run.identifications(report, level, ...))
  
  data = do.call(rbind, lapply(names(run_identifications), function(name) {
    data = reshape(
      run_identifications[[name]],
      direction = 'long', 
      varying = c('complete', 'shared', 'sparse', 'unique'), 
      timevar = 'category', 
      times = c('complete', 'shared', 'sparse', 'unique'), 
      v.names = 'count'
    )
    data$category = factor(
      data$category, 
      levels = rev(c('complete', 'shared', 'sparse', 'unique'))
    )
    data = subset(data, count != 0)
    data$run = as.factor(data$run)
    data$dataset = factor(name, levels = names(run_identifications))
    data
  }))
  
  pl = ggplot(
    data,
    mapping = aes(x = run, y = count)
  ) +
    geom_bar(
      aes(fill = category),
      stat = 'identity',
      alpha = 0.75, color = NA
    ) +
    geom_text(
      aes(
        label = count, color = category, 
        hjust = ifelse(
          count < 20,
          c(
            'complete' = 'middle', 'shared' = 'right', 
            'sparse' = 'left', 'unique' = 'right'
          )[as.character(category)],
          'middle'
        ),
        vjust = c(
          'complete' = 'center', 'shared' = 'center', 
          'sparse' = 'center', 'unique' = 'bottom'
        )[as.character(category)]
      ),
      position = position_stack(vjust = 0.5),
      angle = 0, 
      size = 3
    ) +
    scale_fill_manual(values = c(
      unique = '#ff3355', shared = '#339dff', 
      sparse = '#ffa447', complete = '#65c3ba'
    )) +
    scale_colour_manual(values = c(
      unique = '#7f192a', shared = '#194e7f', 
      sparse = '#997e2a', complete = '#3c756f'
    )) +
    scale_x_discrete(
      name = 'Run'
    ) +
    scale_y_continuous(
      name = list(
        'precursor' = '# Precursors',
        'glycopeptide' = '# Glycopeptides',
        'siteglycan' = '# Site-specific Glycans',
        'glycosite' = '# Protein Glycosites'
      )[[level]],
      expand = expand_scale(mult = c(0, 0.075))
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
    list(plot = pl, data = run_identifications)
  }
  else {
    pl
  }
}


if (exists('reports') && exists('output_dir')) {
  run_identifications = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite', 'peptide', 'protein')
    
    run_identifications = lapply(levels, function(level) {
      res = plot.compare.run.identifications(reports, level, return_data = TRUE)
      
      res$plot = res$plot + theme(strip.text = element_blank())
      
      ggsave(
        paste0(output_dir, '/', 'bar_compare_run_identification_', level, '.svg'), 
        res$plot, 
        width = 12, height = 6, unit = 'cm'
      )
      
      res$data
    })
    names(run_identifications) = levels
    run_identifications
  })
}



plot.compare.cumulative.identifications = function(reports, level, ..., return_data = FALSE) {
  cumulative_identifications = lapply(reports, function(report) {
    get.cumulative.identifications(report, level, ...)
  })
  
  data = do.call(rbind, lapply(names(cumulative_identifications), function(name) {
    data = reshape(
      cumulative_identifications[[name]], 
      direction = 'long', 
      varying = c('full', 'sparse'), 
      timevar = 'category', 
      times = c('full', 'sparse'), 
      v.names = 'count'
    )
    data$dataset = factor(name, levels = names(cumulative_identifications))
    data
  }))
  
  pl = ggplot(
    data,
    mapping = aes(x = run, y = count)
  ) +
    geom_bar(
      aes(fill = category),
      stat = 'identity', 
      position = position_dodge(),
      alpha = 0.75, color = NA
    ) +
    geom_text(
      aes(label = count, color = category),
      position = position_dodge(width = 0.9),
      angle = 90, hjust = 'right',
      size = 3
    ) +
    scale_fill_manual(values = c(
      sparse = '#339dff', full = '#65c3ba'
    )) +
    scale_colour_manual(values = c(
      sparse = '#194e7f', full = '#3c756f'
    )) +
    scale_x_continuous(
      name = 'Cumulative Runs'
    ) +
    scale_y_continuous(
      name = list(
        'precursor' = '# Precursors',
        'glycopeptide' = '# Glycopeptides',
        'siteglycan' = '# Site-specific Glycans',
        'glycosite' = '# Protein Glycosites'
      )[[level]],
      expand = expand_scale(mult = c(0, 0.075))
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
    list(plot = pl, data = cumulative_identifications)
  }
  else {
    pl
  }
}


if (exists('reports') && exists('output_dir')) {
  cumulative_identifications = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite', 'peptide', 'protein')
    
    cumulative_identifications = lapply(levels, function(level) {
      res = plot.compare.cumulative.identifications(reports, level, return_data = TRUE)
      
      res$plot = res$plot + theme(strip.text = element_blank())
      
      ggsave(
        paste0(output_dir, '/', 'bar_compare_cumulative_identification_', level, '.svg'), 
        res$plot, 
        width = 12, height = 6, unit = 'cm'
      )
      
      res$data
    })
    names(cumulative_identifications) = levels
    cumulative_identifications
  })
}


plot.compare.venn.identification = function(reports, level, filename, 
                                    min_frequency = 0.6,
                                    fill = c('#339dff', '#ff3355', '#65c3ba'), 
                                    ...) {
  ident_names = lapply(reports, function(report) {
    identifications = get.identifications(report, level = level)
    
    ident_names = lapply(identifications, function(report) {
      apply(report, 1, paste, collapse = ' ')
    })
    run_counts = table(unlist(ident_names))
    names(run_counts)[run_counts >= length(ident_names) * min_frequency]
  })
  
  VennDiagram::venn.diagram(
    ident_names,
    fill = array(fill, dim = length(ident_names)), 
    alpha = rep(0.5, length(ident_names)),
    # cex = 0, cat.cex = 0,
    fontfamily = "sans", cat.fontfamily = "sans",
    col = 'white', imagetype = 'svg',
    filename = filename, 
    ext.text = FALSE,
    ...
  )
}


if (exists('reports') && exists('output_dir')) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  invisible(local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite', 'peptide', 'protein')
    lapply(levels, function(level) {
      plot.compare.venn.identification(
        rev(reports), level = level,
        filename = paste0(output_dir, '/', 'venn_compare_identification_', level, '.svg'),
        height = 2.5, width = 2.5, unit = 'cm', margin = 0.1
      )
    })
  }))
}

