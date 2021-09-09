arrange.pglyco.report = function(report) {
  
  modified_sequence = sapply(1:nrow(report), function(i) {
    if (is.na(report$Mod[i])) {
      return(report$Peptide[i])
    }
    
    aa = strsplit(report$Peptide[i], '')[[1]]
    
    sapply(strsplit(report$Mod[i], ';')[[1]], function(t) {
      if (t == '') {
        return()
      }
      s = strsplit(t, ',')[[1]]
      pos = as.integer(s[1])
      name = s[2]
      if (name == 'Acetyl[ProteinN-term]') {
        mass = 43
        aa[pos] <<- paste0('n[', mass, ']', aa[pos])
      }
      else {
        if (name == 'Carbamidomethyl[C]') {
          mass = 160
        }
        else if (name == 'Oxidation[M]') {
          mass = 147
        }
        else {
          stop()
        }
        aa[pos] <<- paste0(aa[pos], '[', mass, ']')
      }
    })
    
    paste0(aa, collapse = '')
  })
  
  glycan_composition = sapply(1:nrow(report), function(i) {
    mc = stringr::str_extract_all(report$PlausibleStruct[i], '[A-Za-z0-9]+')[[1]]
    count = table(mc)
    paste0(sapply(names(count)[order(names(count))], function(x) {
      paste0(x, '(', count[x], ')')
    }), collapse = '')
  })
  
  report$Sequence = report$Peptide
  report$ProteinName = sub('^[0-9]+/', '', report$Proteins)
  report$ProteinGlycoSite = sub('^[0-9]+/', '', report$ProSite)
  
  data.frame(
    filename = report$RawName,
    scan = report$Scan,
    RT = report$RT,
    Sequence = report$Peptide,
    FullPeptideName = modified_sequence,
    GlycanStruct = report$PlausibleStruct,
    GlycanComposition = glycan_composition,
    GlycanSite = report$GlySite,
    Charge = report$Charge,
    mz = report$PrecursorMZ,
    Intensity = if ('MonoArea' %in% colnames(report)) report$MonoArea else NA,
    ProteinName = sub('^[0-9]+/', '', report$Proteins),
    ProteinGlycoSite = sub('^[0-9]+/', '', report$ProSite),
    stringsAsFactors = FALSE
  )
}