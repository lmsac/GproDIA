library(ggplot2)


get.cumulative.identifications = function(report, level, run_groups = NULL, ...) {
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
  
  cumulative_numbers = do.call(rbind, lapply(1:length(run_names), function(j) {
    ident_names = ident_names[run_names[[j]]]
    
    cumulative_ident_names = Reduce(
      function(cum, ident_name) {
        list(
          sparse = union(cum$sparse, ident_name),
          full = intersect(cum$full, ident_name)
        )
      },
      ident_names,
      init = list(
        sparse = ident_names[[1]],
        full = ident_names[[1]]
      ),
      accumulate = TRUE
    )[-1]
    
    cumulative_numbers = do.call(rbind, lapply(1:length(cumulative_ident_names), function(run) {
      data.frame(
        group = names(run_names)[j],
        run = run,
        full = length(cumulative_ident_names[[run]]$full),
        sparse = length(cumulative_ident_names[[run]]$sparse),
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  cumulative_numbers
}


plot.cumulative.identifications = function(report, level, ..., return_data = FALSE) {
  cumulative_identifications = get.cumulative.identifications(report, level, ...)
  
  data = reshape(
    cumulative_identifications, 
    direction = 'long', 
    varying = c('full', 'sparse'), 
    timevar = 'category', 
    times = c('full', 'sparse'), 
    v.names = 'count'
  )
  
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
    list(plot = pl, data = cumulative_identifications)
  }
  else {
    pl
  }
}



if (exists('report') && exists('output_dir')) {
  cumulative_identifications = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite', 'peptide', 'protein')
    
    cumulative_identifications = lapply(levels, function(level) {
      res = plot.cumulative.identifications(
        report, level, 
        run_groups = run_groups,
        return_data = TRUE
      )
      
      ggsave(
        paste0(output_dir, '/', 'bar_cumulative_identification_', level, '.svg'), 
        res$plot, 
        width = 6, height = 6, unit = 'cm'
      )
      
      res$data
    })
    names(cumulative_identifications) = levels
    cumulative_identifications
  })
}