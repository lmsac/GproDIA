get.compare.glycoform.per.glycosite = function(reports, shared_ratio = 0) {
  identifications = lapply(reports, function(report) {
    identifications = do.call(rbind, get.identifications(report, level = 'siteglycan'))
    identifications$count = 1
    identifications = aggregate(count ~ ., identifications, sum)
    identifications = subset(identifications, count > length(unique(report$filename)) * shared_ratio)
    identifications
  })
  
  cmp = merge(
    aggregate(GlycanComposition ~ ProteinName + ProteinGlycoSite, identifications[[1]], length), 
    aggregate(GlycanComposition ~ ProteinName + ProteinGlycoSite, identifications[[2]], length), 
    by = c('ProteinName', 'ProteinGlycoSite'), 
    suffixes = paste0('_', names(identifications)[1:2]),
    all = TRUE
  )
  cmp$shared = sapply(1:nrow(cmp), function(i) {
    length(intersect(
      subset(
        identifications[[1]], 
        ProteinName == cmp$ProteinName[i] & 
          ProteinGlycoSite == cmp$ProteinGlycoSite[i]
      )$GlycanComposition, 
      subset(
        identifications[[2]], 
        ProteinName == cmp$ProteinName[i] & 
          ProteinGlycoSite == cmp$ProteinGlycoSite[i]
      )$GlycanComposition
    ))
  })

  cmp = cmp[cmp$shared > 0, ]
  cmp = cmp[order(
    cmp[[paste0('GlycanComposition_', names(identifications)[1])]] - 
      cmp[[paste0('GlycanComposition_', names(identifications)[2])]], 
    -cmp[[paste0('GlycanComposition_', names(identifications)[2])]], 
    -cmp[[paste0('GlycanComposition_', names(identifications)[1])]]
  ), ]
  
  cmp
}


plot.compare.glycoform.per.glycosite = function(reports, return_data = FALSE) {
  cmp = get.compare.glycoform.per.glycosite(reports)
  datasetNames = sub(
    '^GlycanComposition_', '', 
    grep('^GlycanComposition_', colnames(cmp), value = TRUE)
  )
  
  df = rbind(
    data.frame(
      glycosite = paste(cmp$ProteinName, cmp$ProteinGlycoSite),
      count = -cmp$shared + cmp[[paste0('GlycanComposition_', datasetNames[2])]],
      category = datasetNames[2],
      stringsAsFactors = FALSE
    ),
    data.frame(
      glycosite = paste(cmp$ProteinName, cmp$ProteinGlycoSite),
      count = cmp$shared - cmp[[paste0('GlycanComposition_', datasetNames[1])]],
      category = datasetNames[1],
      stringsAsFactors = FALSE
    ),
    data.frame(
      glycosite = paste(cmp$ProteinName, cmp$ProteinGlycoSite),
      count = cmp$shared,
      category = 'shared',
      stringsAsFactors = FALSE
    ),
    data.frame(
      glycosite = paste(cmp$ProteinName, cmp$ProteinGlycoSite),
      count = -cmp$shared,
      category = 'shared_',
      stringsAsFactors = FALSE
    )
  )
  df$category = factor(df$category, levels = unique(df$category))
  df$glycosite = factor(df$glycosite, levels = unique(df$glycosite))
  
  
  pl = ggplot(
    df, 
    aes(x = glycosite, y = count)
  ) +
    geom_bar(
      aes(fill = category),
      stat = 'identity', 
      alpha = 0.6, color = NA
    ) +
    geom_hline(yintercept = 0, color = 'black') +
    scale_fill_manual(values = c( 
      'shared' = '#65c3ba', 'shared_' = '#65c3ba',
      local({
        x = c('#ff3355', '#339dff')
        names(x) = datasetNames[1:2]
        x
      })
    )) +
    scale_y_continuous(
      name = '# Glycoforms'
    ) +
    scale_x_discrete(
      name = 'Protein Glycosite'
    ) +
    theme(
      axis.line.x = element_blank(), 
      axis.line.y = element_line(), 
      panel.background = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(color = 'black'),
      axis.title.y = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = element_blank(),
      legend.position = 'none'
    )
  
  if (return_data) {
    list(plot = pl, data = cmp)
  }
  else {
    pl
  }
}



if (exists('reports') && exists('output_dir')) {
  glycoforms_per_glycosite = local({
    res = plot.compare.glycoform.per.glycosite(reports, return_data = TRUE)
    
    ggsave(
      paste0(output_dir, '/', 'bar_glycoforms_per_glycosite.svg'), 
      res$plot, 
      width = 14, height = 6, unit = 'cm'
    )
    
    res$data
  })
}



plot.compare.glycoform.type.percentages = function(reports, return_data = FALSE) {
  identifications = lapply(reports, function(report) {
    unique(do.call(rbind, get.identifications(report, level = 'siteglycan')))
  })
  
  df = do.call(rbind, lapply(names(identifications), function(dataset) {
    count = sapply(c('^[^FA]+$', '^[^A]*F[^A]*$', '^[^F]*A[^F]*$', 'A.*F'), function(x) 
      length(grep(x, identifications[[dataset]]$GlycanComposition))
    )
    data.frame(
      dataset = dataset,
      type = c('HexNAc+Hex', 'Fuc', 'NeuAc', 'Fuc+NeuAc'), 
      count = sapply(c('^[^FA]+$', '^[^A]*F[^A]*$', '^[^F]*A[^F]*$', 'A.*F'), function(x) 
        length(grep(x, identifications[[dataset]]$GlycanComposition))
      ),
      count = count,
      perenctage = count / sum(count) * 100
    )
  }))
  
  pl = ggplot(
    df,
    aes(x = 1, y = count, fill = type)
  ) + 
    geom_bar(
      stat = "identity", 
      width = 1, colour = NA, alpha = 0.95
    ) + 
    geom_text(
      aes(
        label = sprintf('%s\n%d\n%.0f%%', type, count, perenctage)
      ),
      position = position_stack(vjust = 0.5),
      size = 3
    ) +
    scale_fill_manual(values = c(
      'HexNAc+Hex' = '#c9e3e3', 'Fuc' = '#ffbdbd', 
      'NeuAc' = '#c9c9ff', 'Fuc+NeuAc' = '#ffe8ce'
    )) + 
    coord_polar(theta = "y") +
    theme_void() +
    theme(legend.position = 'none')
  
  if (length(unique(df$dataset)) > 1) {
    pl = pl + 
      facet_grid(cols = vars(dataset)) +
      theme(
        strip.background = element_blank()
      )
  }
  
  if (return_data) {
    list(plot = pl, data = df)
  }
  else {
    pl
  }
}



if (exists('reports') && exists('output_dir')) {
  glycoform_type_percentages = local({
    res = plot.compare.glycoform.type.percentages(reports, return_data = TRUE)
    
    ggsave(
      paste0(output_dir, '/', 'pie_glycoform_types.svg'), 
      res$plot, 
      width = 8, height = 4, unit = 'cm'
    )
    
    res$data
  })
}
