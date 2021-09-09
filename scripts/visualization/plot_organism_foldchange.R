calculate.fold.change = function(report_matrix, run_groups, intensity_column_suffix = '.mzML') {
  report = report_matrix
  
  if ('decoy' %in% colnames(report)) {
    report = subset(report, decoy == 0)
  }
  
  run_names = names(report)[endsWith(names(report), intensity_column_suffix)]
  mean_intensity = sapply(1:length(run_groups), function(i) {
    index = grep(run_groups[[i]], run_names)
    rowMeans(report[, run_names[index]], na.rm = TRUE)
  })
  colnames(mean_intensity) = names(run_groups)
  
  result = apply(combn(length(run_groups), 2), 2, function(x) {
    mean_intensity[, x[2]] / mean_intensity[, x[1]]
  })
  colnames(result) = apply(combn(length(run_groups), 2), 2, function(x) {
    paste0(
      colnames(mean_intensity)[c(x[2], x[1])], 
      collapse = '/'
    )
  })
  
  result = cbind(
    report[, !colnames(report) %in% run_names],
    mean_intensity,
    result
  )
  result
}


get.organism.fold.change = function(report_matrix, organisms) {
  report = report_matrix
  
  organism = rep(NA, nrow(report))
  sapply(organisms, function(org) {
    organism[grep(org, report$ProteinName, ignore.case = TRUE)] <<- org
  })
  report$organism = organism
  
  report = reshape2::melt(
    report,
    measure.vars = grep('/', colnames(report), value = TRUE),
    variable.name = 'group',
    value.name = 'fold_change'
  )
  report
}


plot.organism.fold.change = function(report_matrix, level, run_groups, 
                                     organisms = c('Human', 'Yeast'),
                                     ...,
                                     expected_fold_changes,
                                     return_data = FALSE) {
  fold_changes = calculate.fold.change(
    report_matrix,
    run_groups = run_groups,
    ...
  )
  fold_changes = get.organism.fold.change(
    fold_changes,
    organisms = organisms
  )
  fold_changes = subset(fold_changes, grepl('/S10$', group))
  
  data = fold_changes
  
  pl = ggplot(
    data,
    aes(x = organism, y = fold_change, fill = organism, color = organism)
  ) + 
    stat_summary(
      fun.data = function(x) {
        r = c(
          quantile(x, probs = 0.25) - 1.5 * IQR(x),
          quantile(x, probs = c(0.25, 0.5, 0.75)),
          quantile(x, probs = 0.75) + 1.5 * IQR(x)
        )
        names(r) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
        r
      },
      geom = 'boxplot', position = position_dodge(width = 1), alpha = 0.5
    ) +
    # stat_summary(
    #   fun.y = mean,
    #   geom = 'point', position = position_dodge(width = 1), shape = 18
    # ) +
    stat_summary(
      fun.data = function(x) {
        m = 2 ^ median(x)
        data.frame(
          label = sprintf(ifelse(abs(m) >= 1, '%.2f', '%.3f'), m), 
          y = as.numeric(quantile(x, 0.8))
        )
      },
      geom = 'text', position = position_dodge(width = 1),
      angle = 90, hjust = 0, vjust = 1
    )
  if (!is.null(expected_fold_changes)) {
    expected_fold_changes = local({
      fold_changes = data.frame(expected_fold_changes, check.names = FALSE)
      fold_changes$organism = rownames(fold_changes)
      fold_changes = reshape2::melt(
        fold_changes,
        measure.vars = grep('/', colnames(fold_changes), value = TRUE),
        variable.name = 'group',
        value.name = 'fold_change'
      )
      fold_changes
    })
    
    pl = pl + geom_crossbar(
      data = expected_fold_changes,
      mapping = aes(
        y = fold_change, 
        ymin = fold_change, 
        ymax = fold_change, 
        color = organism
      ),
      width = 1, size = 0.25, linetype = 'dotted',
      position = position_dodge(width = 1)
    )
  } 
  pl = pl +
    scale_fill_manual(values = c('#339dff', '#65c3ba', '#ff3355')) +
    scale_colour_manual(values = c('#194e7f', '#3c756f', '#7f192a')) +
    scale_y_continuous(
      name = paste(
        list(
          'precursor' = 'Precursor',
          'glycopeptide' = 'Glycopeptide',
          'siteglycan' = 'Site-Glycan',
          'glycosite' = 'Glycosite'
        )[[level]], 
        'FC'
      ),
      trans = 'log2'
    ) +
    facet_grid(
      cols = vars(group)
    ) +
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_blank(), 
      panel.background = element_blank(),
      axis.title.y = element_text(color = 'black'),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = 'black'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'none'
    )
  
  if (return_data) {
    list(plot = pl, data = fold_changes)
  }
  else {
    pl
  }
}

organisms = c('Human', 'Yeast')

run_groups = list(
  'S10' = 'serum_yeast.*_1_1',
  'S12' = 'serum_yeast.*_10_12',
  'S15' = 'serum_yeast.*_10_15'
)

expected_fold_changes = list(
  'S12/S10' = c(
    'Human' = (1 / (1 + 1.2)) / (1 / (1 + 1)), 
    'Yeast' = (1.2 / (1 + 1.2)) / (1 / (1 + 1))
  ),
  'S15/S10' = c(
    'Human' = (1 / (1 + 1.5)) / (1 / (1 + 1)), 
    'Yeast' = (1.5 / (1 + 1.5)) / (1 / (1 + 1))
  )
)    

if (exists('report_matrices') && exists('output_dir')) {
  fold_changes = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite')
    
    fold_changes = lapply(levels, function(level) {
      res = plot.organism.fold.change(
        report_matrices[[level]], 
        level = level, 
        run_groups = run_groups, 
        organisms = organisms,
        expected_fold_changes = {
          if (exists('expected_fold_changes')) expected_fold_changes
          else NULL
        },
        return_data = TRUE
      )
      
      #res$plot = res$plot +
      #  coord_cartesian(ylim = c(0.2, 5))
        
      ggsave(
        paste0(output_dir, '/', 'bar_organism_fold_change_', level, '.svg'), 
        res$plot, 
        width = 5, height = 5, unit = 'cm'
      )
      
      res$data
    })
    names(fold_changes) = levels
    fold_changes
  })
}
