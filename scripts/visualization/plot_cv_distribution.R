library(ggplot2)


plot.cv.distribution = function(report_matrix, level, run_groups, ..., return_data = FALSE) {
  cv = calculate.cv(
    report_matrix,
    run_groups = run_groups,
    ...
  )
  
  data = reshape2::melt(
    cv,
    measure.vars = names(run_groups),
    variable.name = 'group',
    value.name = 'cv'
  )
  
  pl = ggplot(
    data, 
    aes(x = cv, color = group, fill = group)
  ) +
    geom_density(alpha = 0.25) +
    geom_vline(
      data = aggregate(cv ~ group, data, median),
      aes(xintercept = cv, color = group),
      linetype = 'dashed'
    ) +
    geom_text(
      data = aggregate(cv ~ group, data, median),
      aes(
        label = sprintf(ifelse(abs(cv) >= 0.01, '%.1f%%', '%.2f%%'), cv * 100), 
        color = group
      ),
      y = Inf, hjust = -1, vjust = 2
    ) +
    
    scale_fill_manual(values = c('#339dff', '#65c3ba', '#ff3355')) +
    scale_colour_manual(values = c('#194e7f', '#3c756f', '#7f192a')) +
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
      strip.text.y = if (length(unique(data$group)) == 1) element_blank()
        else element_text(color = 'black'),
      legend.position = 'none'
    )
  
  if (length(unique(data$group)) > 1) {
    pl = pl + 
      facet_grid(rows = vars(group)) +
      theme(
        strip.background = element_blank()
      )
  }
  
  if (return_data) {
    list(plot = pl, data = cv)
  }
  else {
    pl
  }
}



if (exists('report_matrices') && exists('output_dir') && exists('run_groups')) {
  cv = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite')
    
    cv = lapply(levels, function(level) {
      res = plot.cv.distribution(
        report_matrices[[level]], 
        level = level, 
        run_groups = run_groups, 
        return_data = TRUE
      )
      
      ggsave(
        paste0(output_dir, '/', 'density_cv_', level, '.svg'), 
        res$plot, 
        width = 5, height = 5, unit = 'cm'
      )
      
      res$data
    })
    names(cv) = levels
    cv
  })
  
  # local({
  #   library(openxlsx)
  #   path = file.path(output_dir, paste0(sub('\\..*$', '', basename(report_file)), '.cv.xlsx'))
  #   wb = createWorkbook()
  #   lapply(names(report_matrices), function(level) {
  #     if (level == 'glycopeptide') { return() }
  #     sheet = list(
  #       'precursor' = 'Precursors',
  #       'glycopeptide' = 'Glycopeptides',
  #       'siteglycan' = 'Site-Glycans',
  #       'glycosite' = 'Glycosites'
  #     )[[level]]
  #     addWorksheet(wb, sheet)
  #     writeData(wb, sheet, subset(cv[[level]], rowSums(!is.na(cv[[level]][names(run_groups)])) > 0))
  #     saveWorkbook(wb, file = path, overwrite = TRUE)
  #   })
  #   message(paste('output:', path))
  # })
}
