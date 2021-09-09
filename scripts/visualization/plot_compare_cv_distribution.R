library(ggplot2)


plot.compare.cv.distribution = function(report_matrices, level, run_groups, ..., 
                                        fill = c('#65c3ba', '#ff3355', '#339dff'),
                                        color = c('#3c756f', '#7f192a', '#194e7f'),
                                        return_data = FALSE) {
  cv = lapply(report_matrices, function(report_matrix) {
    calculate.cv(
      report_matrix,
      run_groups = run_groups,
      ...
    )
  })
  data = do.call(rbind, lapply(names(cv), function(name) {
    data = reshape2::melt(
      cv[[name]],
      measure.vars = names(run_groups),
      variable.name = 'group',
      value.name = 'cv'
    )
    data$dataset = factor(name, levels = names(cv))
    data
  }))
  
  pl = ggplot(
    data, 
    aes(x = cv, color = dataset, fill = dataset)
  ) +
    geom_density(alpha = 0.25) +
    geom_vline(
      data = aggregate(cv ~ dataset + group, data, median),
      aes(xintercept = cv, color = dataset),
      linetype = 'dashed'
    ) +
    geom_text(
      data = aggregate(cv ~ dataset + group, data, median),
      aes(
        label = sprintf(ifelse(abs(cv) >= 0.01, '%.1f%%', '%.2f%%'), cv * 100), 
        color = dataset
      ),
      y = Inf, hjust = -0.25, vjust = 2
    ) +
    
    scale_fill_manual(values = fill) +
    scale_colour_manual(values = color) +
    scale_x_continuous(
      name = paste(
        list(
          'precursor' = 'Precursor',
          'glycopeptide' = 'Glycopeptide',
          'siteglycan' = 'Site-Glycan',
          'glycosite' = 'Glycosite'
        )[[level]], 
        'CV'
      )
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    facet_grid(rows = vars(dataset), cols = vars(group)) +
    theme(
      axis.line.x = element_line(),
      axis.line.y = element_blank(),
      panel.background = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(color = 'black'),
      strip.background = element_blank(),
      strip.text.y = {
        if (length(unique(data$dataset)) == 1) element_blank() 
        else element_text(color = 'black')
      },
      strip.text.x = {
        if (length(unique(data$group)) == 1) element_blank()
        else element_text(color = 'black')
      },
      legend.position = 'none'
    )
  
  if (return_data) {
    list(plot = pl, data = cv)
  }
  else {
    pl
  }
}


if (exists('report_matrices_list') && exists('output_dir')) {
  cv = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite')
    
    cv = lapply(levels, function(level) {
      res = plot.compare.cv.distribution(
        lapply(report_matrices_list, function(report_matrices) report_matrices[[level]]), 
        level = level, 
        run_groups = run_groups,
        return_data = TRUE
      )
      
      res$plot = res$plot + theme(strip.text.y = element_blank())
      
      ggsave(
        paste0(output_dir, '/', 'density_compare_cv_', level, '.svg'), 
        res$plot, 
        width = 5, height = 5, unit = 'cm'
      )
      
      res$data
    })
    names(cv) = levels
    cv
  })
  
  # cv_comparison_data = lapply(cv, function(report_list) {
  #   get.comparison.report(report_list, run_groups = run_groups)
  # })
  # 
  # local({
  #   library(openxlsx)
  #   path = file.path(output_dir, paste0(sub('\\..*$', '', basename(report_files[[1]])), '_compare_cv.xlsx'))
  #   wb = createWorkbook()
  #   lapply(names(cv), function(level) {
  #     if (level == 'glycopeptide') { return() }
  #     sheet = list(
  #       'precursor' = 'Precursors',
  #       'glycopeptide' = 'Glycopeptides',
  #       'siteglycan' = 'Site-Glycans',
  #       'glycosite' = 'Glycosites'
  #     )[[level]]
  #     addWorksheet(wb, sheet)
  #     writeData(wb, sheet, cv_comparison_data[[level]])
  #     saveWorkbook(wb, file = path, overwrite = TRUE)
  #   })
  #   message(paste('output:', path))
  # })
}
