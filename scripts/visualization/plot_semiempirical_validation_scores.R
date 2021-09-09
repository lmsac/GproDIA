library(ggplot2)

plot.intensity.similarities = function(scores) {
  pl = ggplot(
    scores, 
    aes(x = dataset, y = intensity_similarity)
  ) + 
    geom_violin(fill = '#339dff') +
    stat_summary(
      fun.data = function(x) {
        r = quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
        names(r) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
        r
      },
      geom = 'boxplot',
      fill = 'white',
      mapping = aes(width = 0.06)
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          label = sprintf('%.3f', median(x)),
          y = 1.05
        )
      },
      geom = 'text'
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          label = sprintf('italic(n)==%d', length(x)),
          y = 0
        )
      },
      geom = 'text',
      hjust = 0, vjust = 1.5,
      angle = 90, parse = TRUE
    ) +
    scale_y_continuous(
      name = 'Dot Product',
      limits = c(0, 1.05),
      breaks = seq(0, 1, 0.2)
    ) +
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_blank(), 
      panel.background = element_blank(),
      axis.title.y = element_text(color = 'black'),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'none'
    )
  
  pl
}

plot.rt.error = function(scores, rt_relative_error = FALSE) {
  pl = ggplot(
    scores, 
    aes(
      x = dataset, 
      y = if (rt_relative_error) { 
        delta_rt * 100 / (
          if ('rt_empirical' %in% columns(scores)) {
            rt_empirical 
          }
          else { rt1 }
        )
      }
      else { delta_rt / 60 }
    )
  ) + 
    stat_summary(
      fun.data = function(x) {
        r = quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
        names(r) = c('ymin', 'lower', 'middle', 'upper', 'ymax')
        r
      },
      geom = 'boxplot', position = position_dodge(), alpha = 0.5,
      color = '#194e7f', fill = '#339dff'
    ) +
    geom_hline(
      aes(yintercept = 0),
      linetype = 'dashed'
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          label = sprintf('IQR=%.1f', IQR(x)),
          y = 5
        )
      },
      geom = 'text',
      hjust = 1, vjust = 1, #nudge_x = 0.05, nudge_y = 0.05,
      angle = 90
    ) +
    stat_summary(
      fun.data = function(x) {
        data.frame(
          label = sprintf('italic(R)[plain("95%%")]*plain("=%.1f")', diff(quantile(x, c(0.025, 0.975)))),
          y = -5
        )
      },
      geom = 'text',
      hjust = 0, vjust = 0, #nudge_x = -0.05, nudge_y = 0.05,
      angle = 90, parse = TRUE
    ) +
    # coord_cartesian(ylim = c(-5, 5)) +
    scale_y_continuous(
      name = if (rt_relative_error) { 'Relative Error (%)' }
      else { 'Delta RT (min)' }
    ) + 
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_blank(), 
      panel.background = element_blank(),
      axis.title.y = element_text(color = 'black'),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'none'
    )
  
  pl
}


plot.rt.correlations = function(scores) {
  get_density = function(x, y, ...) {
    dens = MASS::kde2d(x, y, ...)
    ix = findInterval(x, dens$x)
    iy = findInterval(y, dens$y)
    ii = cbind(ix, iy)
    return(dens$z[ii])
  }
  
  
  rt = do.call(rbind, lapply(unique(scores$dataset), function(x) {
    data = subset(scores, dataset == x)
    rt1 = (
      if ('rt_empirical' %in% colnames(scores)) { data$rt_empirical }
      else { data$rt1 }
    ) / 60
    rt2 = (
      if ('rt_empirical' %in% colnames(scores)) { data$rt_semiempirical }
      else { data$rt2 } 
    ) / 60
    df = data.frame(
      dataset = x,
      rt1 = rt1, rt2 = rt2,
      density = get_density(rt1, rt2, n = 100),
      stringsAsFactors = FALSE
    )
    df = df[order(df$density), ]
    df
  }))
  
  rt_correl = do.call(rbind, lapply(unique(scores$dataset), function(x) {
    data = subset(scores, dataset == x)
    rt1 = (
      if ('rt_empirical' %in% colnames(scores)) { data$rt_empirical }
      else { data$rt1 }
    ) / 60
    rt2 = (
      if ('rt_empirical' %in% colnames(scores)) { data$rt_semiempirical }
      else { data$rt2 } 
    ) / 60
    df = data.frame(
      dataset = x,
      correl = cor(rt1, rt2),
      count = length(rt1),
      stringsAsFactors = FALSE
    )
    df
  }))
  
  pl = ggplot(
    rt, 
    aes(x = rt1, y = rt2, color = density)
  ) +
    geom_point(shape = 16, alpha = 0.75, size = 0.5) +
    scale_color_gradient2(
      high = '#0a1f33', low = '#c1e1ff', mid = '#194e7f',
      space = 'Lab', midpoint = quantile(rt$density, 0.75),
      name = 'Density'
    ) +
    scale_x_continuous(
      name = if ('rt_empirical' %in% colnames(scores)) {
        'Experimental RT (min)'
      } else {
        'RT 1 (min)'
      }
    ) +
    scale_y_continuous(
      name = if ('rt_semiempirical' %in% colnames(scores)) {
        'Predicted RT (min)'
      } else {
        'RT 2 (min)'
      }
    ) +
    
    geom_text(
      data = rt_correl,
      mapping = aes(
        label = sprintf(
          'atop(italic(r)*plain("=%.4f")*plain("\n"), italic(n)*plain("=%d"))', 
          correl, count
        )
      ),
      x = min(rt$rt1), y = max(rt$rt2),
      color = 'black',
      hjust = 0, vjust = 1,
      parse = TRUE
    ) +
    
    #coord_cartesian(xlim = c(10, 50), ylim = c(10, 50)) +
    facet_wrap(vars(dataset)) +
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_line(), 
      panel.background = element_blank(),
      axis.title.y = element_text(color = 'black'),
      axis.title.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = element_text(color = 'black'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'none',
      aspect.ratio = 1
    )
  
  pl
}



cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) >= 2) {
  score_files = cmdArgs[seq(2, length(cmdArgs), 2)]
  names(score_files) = cmdArgs[seq(2, length(cmdArgs), 2) - 1]
}

if (length(cmdArgs) %% 2 == 1) {
  output_dir = tail(cmdArgs, 1)
} else {
  output_dir = 'plots'
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# score_files = list.files(
#   'fissionyeast/library/semiempirical_validation',
#   pattern = 'semiempirical_validation_.*\\.score\\.csv$',
#   full.names = TRUE
# )
# names(score_files) = local({
#   s = sub('^.*semiempirical_validation_(.+)\\.score\\.csv$', '\\1', score_files)
#   ifelse(
#     grepl('^k[0-9]+$', s), substring(s, 2),
#     ifelse(s == 'noknn', 'no KNN', s)
#   )
# })


library(readr)

scores = lapply(score_files, read_csv)

scores = do.call(
  rbind, 
  lapply(names(scores), function(x) {
    data.frame(
      dataset = x,
      scores[[x]],
      stringsAsFactors = FALSE
    )
  })
)
scores$dataset = factor(scores$dataset, levels = unique(scores$dataset))

local({
  pl = plot.intensity.similarities(scores)
  pl = pl + 
    scale_x_discrete(
      name = '# Nearest Neighbors'
    ) +
    theme(
      axis.title.x = element_text(color = 'black')
    )
  
  ggsave(
    file.path(output_dir, 'intensity_similarity.svg'), 
    pl, 
    width = 8, height = 6, unit = 'cm'
  )
})


local({
  pl = plot.rt.error(scores)
  pl = pl + 
    scale_x_discrete(
      name = '# Nearest Neighbors'
    ) +
    coord_cartesian(ylim = c(-5, 5)) +
    theme(
      axis.title.x = element_text(color = 'black')
    )
  
  ggsave(
    file.path(output_dir, 'rt_error.svg'), 
    pl,
    width = 8, height = 6, unit = 'cm'
  )
})


local({
  pl = plot.rt.correlations(scores)
  pl = pl + 
    coord_cartesian(xlim = c(2, 60), ylim = c(2, 60))
  
  ggsave(
    file.path(output_dir, 'rt_correlation.svg'), 
    pl,
    width = 16, height = 12, unit = 'cm'
  )
})

