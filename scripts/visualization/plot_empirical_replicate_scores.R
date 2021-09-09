# score_files = list(
#   'Budding Yeast' = 'yeast/library/yeast/yeast_6h_consensus_rtcalibrated.score.csv',
#   'Fission Yeast' = 'fissionyeast/library/fissionyeast_lab/fissionyeast_lab_consensus_rtcalibrated.score.csv',
#   'Serum' = 'serum_/library/serum_fraction/serum_fraction_consensus_rtcalibrated.score.csv'
# )

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
  pl = plot.intensity.similarities(
    aggregate(
      intensity_similarity ~ 
        dataset + peptideSequence + modification + 
        glycanStruct + glycanSite + precursorCharge, 
      scores,
      min
    )
  )
  
  ggsave(
    file.path(output_dir, 'intensity_similarity.svg'), 
    pl, 
    width = 8, height = 6, unit = 'cm'
  )
})

local({
  pl = plot.rt.error(
    aggregate(
      delta_rt ~ 
        dataset + peptideSequence + modification + 
        glycanStruct + glycanSite + precursorCharge, 
      scores,
      function(x) x[which.max(abs(x))]
    )
  )
  
  ggsave(
    file.path(output_dir, 'rt_error.svg'), 
    pl,
    width = 8, height = 6, unit = 'cm'
  )
})


local({
  pl = plot.rt.correlations(
    merge(
      scores, 
      aggregate(cbind(file1, file2) ~ dataset, scores, head, 1),
      all = FALSE
    )
  )
  pl = pl + 
    coord_cartesian(xlim = c(2, 60), ylim = c(2, 60)) +
    theme(strip.text = element_blank())
  
  ggsave(
    file.path(output_dir, 'rt_correlation.svg'), 
    pl,
    width = 15, height = 5, unit = 'cm'
  )
})