library(ggplot2)


get.run.identifications = function(report, level, run_groups = NULL, ...) {
  identifications = get.identifications(report, level = level, ...)
  ident_names = lapply(identifications, function(report) {
    apply(report, 1, paste, collapse = ' ')
  })
  
  if (!is.null(run_groups)) {
    run_names = lapply(run_groups, function(run_group) {
      grep(run_group, names(ident_names), value = TRUE)
    })
  }
  else {
    run_names = list('All' = names(ident_names))
  }
  
  identification_numbers = do.call(rbind, lapply(1:length(run_names), function(j) {
    ident_names = ident_names[run_names[[j]]]
    
    run_counts = table(unlist(ident_names))
    identification_numbers = do.call(rbind, lapply(1:length(ident_names), function(i) {
      run_count = run_counts[ident_names[[i]]]
      data.frame(
        group = names(run_names)[j],
        run = i,
        filename = names(ident_names)[i],
        complete = sum(run_count == length(ident_names)),
        shared = sum(run_count > length(ident_names) * 0.5 & run_count < length(ident_names)),
        sparse = sum(run_count > 1 & run_count <= length(ident_names) * 0.5),
        unique = sum(run_count == 1),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  identification_numbers
}


plot.run.identifications = function(report, level, ..., return_data = FALSE) {
  run_identifications = get.run.identifications(report, level, ...)
  
  data = reshape(
    run_identifications,
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
  
  if (length(unique(data$group)) > 1) {
    pl = pl + 
      facet_grid(cols = vars(group))
  }
  
  if (return_data) {
    list(plot = pl, data = run_identifications)
  }
  else {
    pl
  }
}



if (exists('report') && exists('output_dir')) {
  run_identifications = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite', 'peptide', 'protein')
    
    run_identifications = lapply(levels, function(level) {
      res = plot.run.identifications(
        report, level, 
        run_groups = run_groups,
        return_data = TRUE
      )
      
      ggsave(
        paste0(output_dir, '/', 'bar_run_identification_', level, '.svg'), 
        res$plot, 
        width = 6, height = 6, unit = 'cm'
      )
      
      res$data
    })
    names(run_identifications) = levels
    run_identifications
  })
}

