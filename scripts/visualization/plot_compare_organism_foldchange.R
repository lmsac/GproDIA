plot.compare.organism.fold.change = function(report_matrices, level, run_groups, 
                                     organisms = c('Human', 'Yeast'),
                                     ...,
                                     expected_fold_changes,
                                     return_data = FALSE) {
  fold_changes = lapply(report_matrices, function(report_matrix) {
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
  })
  
  data = do.call(rbind, lapply(names(fold_changes), function(name) {
    data = fold_changes[[name]]
    data$dataset = factor(name, levels = names(fold_changes))
    data
  }))
  
  pl = ggplot(
    data,
    aes(x = dataset, y = fold_change, fill = organism, color = organism)
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
    
    expected_fold_changes = do.call(rbind, lapply(unique(data$dataset), function(dataset) {
      expected_fold_changes$dataset = dataset
      expected_fold_changes
    }))
    
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
      axis.ticks.x = element_line(),
      axis.text.x = element_text(color = 'black'),
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

   

if (exists('report_matrices_list') && exists('output_dir')) {
  fold_changes = local({
    levels = c('precursor', 'glycopeptide', 'siteglycan', 'glycosite')
    
    fold_changes = lapply(levels, function(level) {
      res = plot.compare.organism.fold.change(
        lapply(report_matrices_list, function(report_matrices) report_matrices[[level]]), 
        level = level, 
        run_groups = run_groups, 
        organisms = organisms,
        expected_fold_changes = {
          if (exists('expected_fold_changes')) expected_fold_changes
          else NULL
        },
        return_data = TRUE
      )
      
      # res$plot = res$plot +
      #   scale_y_continuous(
      #     name = paste(
      #       list(
      #         'precursor' = 'Precursor',
      #         'glycopeptide' = 'Glycopeptide',
      #         'siteglycan' = 'Site-Glycan',
      #         'glycosite' = 'Glycosite'
      #       )[[level]], 
      #       'FC'
      #     ),
      #     trans = 'log2',
      #     breaks = c(0.25, 0.5, 1, 2, 4)
      #   ) + 
      #   coord_cartesian(ylim = c(0.25, 4))
      
      ggsave(
        paste0(output_dir, '/', 'bar_compare_organism_fold_change_', level, '.svg'), 
        res$plot, 
        width = 10, height = 6.25, unit = 'cm'
      )
      
      res$data
    })
    names(fold_changes) = levels
    fold_changes
  })
}
