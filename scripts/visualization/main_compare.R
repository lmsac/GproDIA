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

sourceInCurrentDir(c(
  'report.R', 
  'compare_report.R'
))
sourceInCurrentDir(c(
  'plot_run_identifications.R', 
  'plot_cumulative_identifications.R'
))

cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) >= 2) {
  report_files = cmdArgs[seq(2, length(cmdArgs), 2)]
  names(report_files) = cmdArgs[seq(2, length(cmdArgs), 2) - 1]
}

if (length(cmdArgs) %% 2 == 1) {
  output_dir = tail(cmdArgs, 1)
} else {
  output_dir = 'plots'
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# report_files = rev(list(
#   'DIA 1h' = 'fissionyeast/result/fissionyeast/fissionyeast_aligned.tsv',
#   'DDA 1h' = 'fissionyeast/DDA/fissionyeast_1h/pGlycoDB-GP-FDR-Pro-Quant.txt',
#   'DDA 6h' = 'fissionyeast/DDA/fissionyeast_6h/pGlycoDB-GP-FDR-Pro-Quant.txt'
# ))
# 
# report_files = rev(list(
#   'SSL' = 'fissionyeast/result/fissionyeast/fissionyeast_aligned.tsv',
#   'DDA 1h' = 'fissionyeast/DDA/fissionyeast_1h/pGlycoDB-GP-FDR-Pro-Quant.txt',
#   'LRL' = 'fissionyeast/result/fissionyeast_lab/fissionyeast_lab_aligned.tsv'
# ))
# 
# report_files = list(
#   'SSL' = 'fissionyeast/result/fissionyeast/fissionyeast_aligned.tsv',
#   'LRL' = 'fissionyeast/result/fissionyeast_lab/fissionyeast_lab_aligned.tsv',
#   'EXL' = 'fissionyeast/result/fissionyeast_extend/fissionyeast_extend_aligned.tsv'
# )
# 
# report_files = list(
#   'Both NoGF' = 'fissionyeast/result/entrapment_human/fissionyeast_entrapment_human_uis_nogf_nogp.tsv',
#   'Both GF' = 'fissionyeast/result/entrapment_human/fissionyeast_entrapment_human_uis_nogp.tsv',
#   'Peptide NoGF.' = 'fissionyeast/result/entrapment_peptide/fissionyeast_entrapment_peptide_uis_nogf_nogp.tsv',
#   'Peptide GF.' = 'fissionyeast/result/entrapment_peptide/fissionyeast_entrapment_peptide_uis_nogp.tsv',
#   'Glycan NoGF' = 'fissionyeast/result/entrapment_glycan/fissionyeast_entrapment_glycan_uis_nogf_nogp.tsv',
#   'Glycan GF' = 'fissionyeast/result/entrapment_glycan/fissionyeast_entrapment_glycan_uis_nogp.tsv'
# )




library(readr)

if (!exists('report_files')) {
  stop('No report file specified')
}

reports = lapply(report_files, function(report_file){
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
  
  report
})


run_groups = list('All' = '(.*?)')

report_matrices_list = lapply(reports, function(report) {
  get.report.matrices(report, intensity = 'ms2')
})


sourceInCurrentDir(c(
  'plot_compare_identifications.R',
  'plot_compare_cv_distribution.R'
))



