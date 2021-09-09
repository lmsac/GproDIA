library(ggplot2)

plot.rt.calibration = function(rtcalibration_report) {
  pl = ggplot(
    rtcalibration_report
  ) +
    geom_point(
      aes(
        x = rt_old  / 60, 
        y = rt_new / 60
      ), 
      alpha = 0.25, size = 0.05, color = '#339dff'
    ) +
    geom_point(
      aes(
        x = rt_old  / 60, 
        y = rt_reference / 60
      ), 
      size = 2, color = '#ff3355', 
      shape = 18
    ) +
    
    scale_x_continuous(
      name = 'Raw RT (min)'
    ) +
    scale_y_continuous(
      name = 'Calibrated RT (min)'
    ) +
    facet_wrap(vars(run)) +
    
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
if (length(cmdArgs) > 0) {
  rtcalibration_report_file = cmdArgs[1]
}

if (length(cmdArgs) > 1) {
  output_dir = cmdArgs[2]
} else {
  output_dir = 'plots'
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# rtcalibration_report_file = "fissionyeast/library/fissionyeast_6h/fissionyeast_6h_consensus_rtcalibrated.csv"

library(readr)

if (!exists('rtcalibration_report_file')) {
  stop('No report file specified')
}

rtcalibration_report = read_csv(rtcalibration_report_file)
rtcalibration_report$run = paste('Run ', as.integer(factor(rtcalibration_report$run)))

local({
  pl = plot.rt.calibration(rtcalibration_report)
  ggsave(
    paste0(output_dir, '/', 'scatter_rtcalibrartion', '.svg'), 
    pl, 
    width = 14, height = 6, unit = 'cm'
  )
})

