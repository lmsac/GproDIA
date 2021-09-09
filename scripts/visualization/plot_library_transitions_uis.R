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


plot.transitions.uis = function(identifying_transitions, detecting_transitions) {
  transitions = local({
    transitions = identifying_transitions
    
    transitions$glycoform_composition = sapply(strsplit(transitions$glycoform, '[()]'), function(s) {
      s = s[s != '']
      count = table(s)
      
      paste0(paste0(names(count), '(', count, ')'), collapse = '')
    })
    
    glycoforms = aggregate(annotation ~ glycoform + glycoform_composition, transitions, length)
    colnames(glycoforms)[3] = 'num_transitions'
    glycoforms = glycoforms[order(glycoforms$num_transitions, decreasing = TRUE), ]
    glycoforms$num_transitions = NULL
    
    glycan_composition = glycoforms$glycoform_composition[
      match(transitions$glycan[1], glycoforms$glycoform)
      ]
    glycoforms = subset(
      glycoforms, 
      glycoform_composition != glycan_composition & 
        !duplicated(glycoform_composition)
    )
    glycoforms = rbind(
      data.frame(
        glycoform = transitions$glycan[1], 
        glycoform_composition = glycan_composition,
        stringsAsFactors = FALSE
      ), 
      glycoforms
    )
    glycoforms = head(glycoforms, 3)
    
    transitions = do.call(rbind, lapply(1:nrow(glycoforms), function(i) {
      subset(transitions, glycoform == glycoforms$glycoform[i])
    }))
    
    transitions$shared = local({
      transition_glycoform_count = table(transitions$annotation) 
      transition_glycoform_count = transition_glycoform_count[transitions$annotation]
      
      ifelse(
        transition_glycoform_count == length(unique(transitions$glycoform)), 'common',
        ifelse(transition_glycoform_count == 1, 'unique', 'shared')            
      )
    })
    
    transitions
  })
  
  detecting_transitions = subset(
    detecting_transitions,
    fragment_type != 'Y' | (annotation %in% transitions$annotation[transitions$shared == 'common'])
  )
  identifying_transitions = subset(
    transitions, 
    (product_charge %in% c(2, 3)) & 
      fragment_type == 'Y' &
      !(annotation %in% detecting_transitions$annotation) &
      shared != 'common' &
      sapply(product_mz, function(x) sum(abs(detecting_transitions$product_mz - x) < 5) == 0)
  )
  identifying_transitions$glycoform_composition = factor(
    identifying_transitions$glycoform_composition,
    levels = unique(identifying_transitions$glycoform_composition)
  )
  
  annotation.to.label = function(annotation) {
    ifelse(
      annotation == 'MS2_Precursor_i0', 'Prec',
      local({
        t = gsub('^DECOY_', 'D', annotation)
        t = gsub('\\[[+\\-][0-9.]+\\]$', '', t)
        t = gsub('[()]', '', t)
        t = sub('\\$', '*plain("$")', t)
        t = gsub('([^+\\-])([0-9]+)', '\\1[\\2]', t)
        t = gsub('\\]([A-Za-z])', ']*\\1', t)
        t = sub('^(Y-[^\\^]+)\\^([^\\^]+)', '{\\1}^\\2', t)
        t
      })
    )
  }
  pl = ggplot() +
    geom_linerange(
      data = detecting_transitions,
      mapping = aes(
        x = product_mz, ymin = 0, 
        ymax = product_intensity / max(product_intensity, na.rm = TRUE) * 100,
        color = 'common'
      )
    ) +
    geom_linerange(
      data = identifying_transitions,
      mapping = aes(
        x = product_mz, ymin = 0, 
        ymax = 100,
        color = ifelse(shared == 'common', 'common', glycoform)
      ),
      linetype = 'dashed'
    ) +
    geom_text(
      data = cbind(
        detecting_transitions, 
        glycoform_composition = local({
          count = table(identifying_transitions$glycoform_composition)
          names(count)[which.min(count)]
        })
      ),
      mapping = aes(
        x = product_mz, 
        y = product_intensity / max(product_intensity, na.rm = TRUE) * 100,
        label = annotation.to.label(annotation),
        hjust = ifelse(
          product_intensity >= max(product_intensity, na.rm = TRUE) * 0.55, 
          'right', 'left'
        ),
        vjust = ifelse(
          product_intensity >= max(product_intensity, na.rm = TRUE) * 0.55, 
          'top', 'center'
        )
      ),
      angle = 90, size = 2,
      parse = TRUE
    ) +
    geom_text(
      data = identifying_transitions,
      mapping = aes(
        x = product_mz, 
        y = ifelse(
          !is.na(product_intensity), 
          product_intensity / max(product_intensity, na.rm = TRUE) * 100,
          100
        ),
        label = annotation.to.label(annotation), 
        hjust = 'right',
        vjust = 'top'
      ),
      angle = 90, size = 2,
      parse = TRUE
    ) +
    scale_x_continuous(
      name = parse(text = 'italic(m)/italic(z)'),
      limits = c(750, 1850)
    ) +
    scale_y_continuous(
      name = 'Intensity'
    ) +
    scale_color_manual(values = c(
      '#ffa447', '#ff3355', '#339dff', '#65c3ba'
    )) +
    facet_grid(
      rows = vars(glycoform_composition), 
      labeller = function(x) {
        t = gsub('\\(', '[', x$glycoform_composition)
        t = gsub('\\)', ']', t)
        t = gsub('\\]([A-Za-z])', ']*\\1', t)
        x$glycoform_composition = t
        label_parsed(x)
      }) +
    theme(
      axis.line.y = element_line(), 
      axis.line.x = element_line(), 
      panel.background = element_blank(),
      axis.title.y = element_text(color = 'black'),
      axis.title.x = element_text(color = 'black'),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = 'black'),
      strip.text = element_text(face = 'bold'),
      legend.position = 'none'
    )
  
  pl
}



cmdArgs = commandArgs(trailingOnly = TRUE)
if (length(cmdArgs) >= 5) {
  osw_file = cmdArgs[1]
  modified_sequence = cmdArgs[2]
  glycan_site = as.integer(cmdArgs[3])
  glycan_name = cmdArgs[4]
  precursor_charge = cmdArgs[5]
}

if (length(cmdArgs) >= 6) {
  output_dir = cmdArgs[6]
} else {
  output_dir = 'plots'
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}



# osw_file = 'fissionyeast/result/entrapment_glycan/fissionyeast_entrapment_glycan_uis.osw'
# modified_sequence = 'VYAVJASLNHDK'
# glycan_name = '(N(N(H(H(H(H)))(H(H(H(H(H(H)))))(H(H(H(H(H)))))))))'
# glycan_site = 5
# precursor_charge = 3


if (exists('osw_file') && exists('output_dir')) {
  transitions = local({
    identifying_transitions = read.identifying.transitions(
      osw_file,
      modified_sequence,
      glycan_site,
      glycan_name,
      precursor_charge
    )
    detecting_transitions = read.detecting.transitions(
      osw_file,
      modified_sequence,
      glycan_site,
      glycan_name,
      precursor_charge
    )
    
    pl = plot.transitions.uis(identifying_transitions, detecting_transitions)
    # pl = pl + theme(
    #   strip.text = element_blank()
    # )
    
    ggsave(
      paste0(
        output_dir, '/', 
        sprintf(
          'transitions_uis_%s_%s_%s_%s.svg', 
          modified_sequence,
          glycan_site,
          glycan_name,
          precursor_charge
        )
      ), 
      plot = pl,
      width = 8, height = 8, unit = 'cm'
    )
    
    list(
      identifying_transitions = identifying_transitions,
      detecting_transitions = detecting_transitions
    )
  })
}

