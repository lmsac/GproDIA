sourceInCurrentDir = local({
  getCurrentScript = function() {
    if (Sys.getenv("RSTUDIO") == "1" && rstudioapi::isAvailable(version_needed = NULL)) {
      # RStudio interactive
      return(normalizePath(rstudioapi::getSourceEditorContext()$path))
    }
    
    lastScriptSourced = tail(unlist(lapply(sys.frames(), function(env) env$ofile)), 1)
    if (!is.null(lastScriptSourced)) {
      # 'source'd via R console
      return(normalizePath(lastScriptSourced, winslash = .Platform$file.sep, mustWork = TRUE))
    }
    
    cmdArgs = commandArgs(trailingOnly = FALSE)
    
    # Rscript/R console option
    needle = "--file="
    match = grep(needle, cmdArgs)
    if (length(match) > 0) {
      return(normalizePath(sub(needle, "", cmdArgs[match]), winslash = .Platform$file.sep, mustWork = TRUE)[1])
    }
    
    # R console option
    match = grep("^-f$", cmdArgs)
    if (length(match) > 0) {
      return(normalizePath(dirname(cmdArgs[match + 1]))[1])
    }
    
    NULL
  }
  
  currentScript = getCurrentScript()
  
  sourceInCurrentDir = function(file, ...) {
    invisible(lapply(file.path(dirname(currentScript), file), source, ...))
  }
  
  sourceInCurrentDir
})

sourceInCurrentDir(c('report.R'))


cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) > 0) {
  report_file = cmdArgs[1]
}

if (length(cmdArgs) > 1) {
  output_dir = cmdArgs[2]
} else {
  output_dir = 'report'
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# report_file = 'fissionyeast/result/fissionyeast/fissionyeast_aligned.tsv'

library(readr)

if (!exists('report_file')) {
  stop('No report file specified')
}

report = read_delim(
  report_file, 
  '\t', escape_double = FALSE, trim_ws = TRUE,
  col_types = cols(ProteinGlycoSite = col_character())
)

if (any(c('RawName', 'PlausibleStruct') %in% colnames(report))) {
  if (!exists('arrange.pglyco.report')) {
    sourceInCurrentDir('pglyco_report.R')
  }
  report = arrange.pglyco.report(report)
  report$filename = paste0(report$filename, '.mzML')
}

run_groups = list('All' = '(.*?)')

report_matrices = get.report.matrices(report, intensity = 'ms2')


local({
  library(openxlsx)
  path = file.path(output_dir, paste0(sub('\\..*$', '', basename(report_file)), '.matrix.xlsx'))
  wb = createWorkbook()
  lapply(names(report_matrices), function(level) {
    if (level == 'glycopeptide') { return() }
    sheet = list(
      'precursor' = 'Precursors',
      'glycopeptide' = 'Glycopeptides',
      'siteglycan' = 'Site-Glycans',
      'glycosite' = 'Glycosites'
    )[[level]]
    addWorksheet(wb, sheet)
    writeData(wb, sheet, report_matrices[[level]])
    saveWorkbook(wb, file = path, overwrite = TRUE)
  })
  message(paste('output:', path))
})

# lapply(names(report_matrices), function(level) {
#   path = file.path(output_dir, paste0(sub('\\.tsv$', '', basename(report_file)), '.', level, '_matrix.csv'))
#   write_csv(
#     report_matrices[[level]],
#     path = path,
#     na = ''
#   )
#   message(paste('output:', path))
# })


sourceInCurrentDir(c(
  'plot_run_identifications.R',
  'plot_cumulative_identifications.R',
  'plot_cv_distribution.R'
))
