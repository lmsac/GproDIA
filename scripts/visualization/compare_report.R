library(reshape2)

get.comparison.report = function(report_list, run_groups, 
                                 group.name = 'group', value.name = 'value') {
  data = do.call(rbind, lapply(
    names(report_list), 
    function(label) cbind(
      dataset = label, 
      report_list[[label]]
    )
  ))
  
  if (group.name %in% colnames(data)) {
    data = data[, setdiff(colnames(data), names(run_groups)), drop = FALSE]
  }
  else {
    data = melt(
      data, 
      measure.vars = names(run_groups), 
      variable.name = group.name,
      value.name = value.name,
    )
  }
  
  data = subset(data, !is.na(data[[value.name]]))
  
  data = dcast(
    data,  
    formula = paste0(
      paste(setdiff(colnames(data), c('dataset', group.name, value.name)), collapse = '+'),
      ' ~ dataset + ', group.name
    ), 
    value.var = value.name
  )
  data
}












