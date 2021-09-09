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
  
  sourceInCurrentDir(c('library.R'))
  
  sourceInCurrentDir
})


library(ggplot2)

plot.transitions = function(transitions) {
  pl = ggplot(transitions) +
    geom_linerange(
      aes(x = product_mz, ymin = 0, 
          ymax = product_intensity / max(product_intensity) * 100,
          color = fragment_type)
    ) +
    geom_text(
      aes(
        x = product_mz, y = product_intensity / max(product_intensity) * 100,
        label = local({
          t = gsub('^DECOY_', 'D', annotation)
          t= gsub('\\[[+\\-][0-9.]+\\]$', '', t)
          t= gsub('[()]', '', t)
          t = sub('\\$', '*plain("$")', t)
          t= gsub('([^+\\-])([0-9]+)', '\\1[\\2]', t)
          t= gsub('\\]([A-Za-z])', ']*\\1', t)
          t= sub('^(Y-[^\\^]+)\\^([^\\^]+)', '{\\1}^\\2', t)
        }),
        hjust = ifelse(product_intensity == max(product_intensity), 'right', 'left'),
        vjust = ifelse(product_intensity == max(product_intensity), 'top', 'center')
      ),
      angle = 90, size = 2,
      parse = TRUE
    ) +
    scale_x_continuous(
      name = parse(text = 'italic(m)/italic(z)'),
      limits = c(200, 1600)
    ) +
    scale_y_continuous(
      name = 'Intensity (%)'
    ) +
    scale_color_manual(values = c(
      'b' = '#339dff', 
      'y' = '#339dff', 
      'b-N(1)' = '#65c3ba',
      'y-N(1)' = '#65c3ba',
      'b$' = '#65c3ba',
      'y$' = '#65c3ba',
      'Y' = '#ff3355',
      'DECOY_Y' = '#ffa447'
    )) +
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_line(), 
      panel.background = element_blank(),
      axis.title.y = element_text(color = 'black'),
      axis.title.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black'),
      axis.text.x = element_text(color = 'black'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'none'
    )
  
  pl
}


get.decoy = function(modified_sequence, glycan_site, glycan_name,
                     decoy_peptide = TRUE, decoy_glycan = TRUE) {
  if (decoy_peptide) {
    modified_sequence = stringr::str_extract_all(modified_sequence, '[A-Z](\\[[0-9]+\\])?')[[1]]
    aa_length = length(modified_sequence)
    if (modified_sequence[aa_length] %in% c('K', 'R')) {
      modified_sequence = c(
        rev(modified_sequence[-aa_length]), 
        modified_sequence[aa_length]
      )
      glycan_site = aa_length - glycan_site
    } else {
      modified_sequence = rev(modified_sequence)
      glycan_site = aa_length - glycan_site + 1
    }
    modified_sequence = paste0(modified_sequence, collapse = '')
  }
  if (decoy_glycan) {
    glycan_name = paste0('DECOY_', glycan_name)
  }
  
  list(
    modified_sequence = modified_sequence,
    glycan_site = glycan_site, 
    glycan_name = glycan_name
  )
}


cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) >= 5) {
  pqp_file = cmdArgs[1]
  modified_sequence = cmdArgs[2]
  glycan_site = as.integer(cmdArgs[3])
  glycan_name = cmdArgs[4]
  precursor_charge = cmdArgs[5]
}

if (length(cmdArgs) >= 6) {
  if (tolower(cmdArgs[6]) == 'decoy_peptide') {
    decoy_peptide = TRUE
    decoy_glycan = FALSE
  }
  else if (tolower(cmdArgs[6]) == 'decoy_glycan') {
    decoy_peptide = FALSE
    decoy_glycan = TRUE
  }
  else if (tolower(cmdArgs[6]) == 'decoy_both') {
    decoy_peptide = TRUE
    decoy_glycan = TRUE
  }
  else {
    decoy_peptide = FALSE
    decoy_glycan = FALSE
  }
} 

if (get0('decoy_peptide', ifnotfound = FALSE) | get0('decoy_glycan', ifnotfound = FALSE)) {
  if (length(cmdArgs) >= 7) {
    output_dir = cmdArgs[7]
  } 
  else {
    output_dir = 'plots'
  }
} else {
  if (length(cmdArgs) >= 6) {
    output_dir = cmdArgs[6]
  } 
  else {
    output_dir = 'plots'
  }
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# pqp_file = 'fissionyeast/library/fissionyeast/fissionyeast.PQP'
# modified_sequence = 'VYAVJASLNHDK'
# glycan_name = '(N(N(H(H(H(H)))(H(H(H(H(H(H)))))(H(H(H(H(H)))))))))'
# glycan_site = 5
# precursor_charge = 3


if (exists('decoy_peptide') || exists('decoy_glycan')) {
  local({
    res = get.decoy(
      modified_sequence,
      glycan_site,
      glycan_name,
      decoy_peptide = get0('decoy_peptide', ifnotfound = FALSE),
      decoy_glycan = get0('decoy_glycan', ifnotfound = FALSE)
    )
    
    modified_sequence <<- res$modified_sequence
    glycan_site <<- res$glycan_site 
    glycan_name <<- res$glycan_name
  })
}


if (exists('pqp_file') && exists('output_dir')) {
  transitions = local({
    transitions = read.detecting.transitions(
      pqp_file,
      modified_sequence,
      glycan_site,
      glycan_name,
      precursor_charge
    )
    
    pl = plot.transitions(transitions)
    
    # pl = pl +
    #   theme(
    #     axis.line.y = element_line(), 
    #     axis.title.y = element_blank(),
    #     axis.text.y = element_blank(),
    #     axis.ticks.y = element_blank(),
    #     axis.title.x = element_blank()
    #   )
    
    ggsave(
      paste0(
        output_dir, '/', 
        sprintf(
          'transitions_%s_%s_%s_%s.svg', 
          modified_sequence,
          glycan_site,
          glycan_name,
          precursor_charge
        )
      ), 
      plot = pl,
      width = 6, height = 3, unit = 'cm'
    )
    
    transitions
  })
}

