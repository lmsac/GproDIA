from assay import GlycoAssayBuilder 
from assay.annotation import SpectrumAnnotator
from .pglyco2assay import pGlycoToAssayConverter

import numpy as np
import pandas as pd


def remove_suspicious_glycan_struct(psm_report):
    composition_column = psm_report.columns.values \
        [psm_report.columns.str.contains('^Glycan\\(.*\\)$')][0]
    
    def find_best_struct(x):
        score = x.groupby('PlausibleStruct')['TotalScore'].sum() \
            .sort_values(ascending=False)
        struct = score.index[0]
        return x.loc[x['PlausibleStruct'] == struct]
    
    data = pd.concat(
        (psm_report, pd.DataFrame(
            list(range(0, len(psm_report))), 
            columns=['index'], 
            index=psm_report.index
        )), 
        axis=1
    ).fillna('NA')
    
    data = data \
        .groupby(
            ['Peptide', 'Mod', 'Charge', 'GlySite', composition_column]
        ) \
        .apply(find_best_struct) \
        .reset_index(drop=True)
        
    return psm_report.iloc[data['index'].sort_values()]
        

def remove_suspicious_glycan_site(psm_report):
    def find_best_site(x):
        score = x.groupby('GlySite')['TotalScore'].sum() \
            .sort_values(ascending=False)
        site = score.index[0]
        return x.loc[x['GlySite'] == site]
    
    data = pd.concat(
        (psm_report, pd.DataFrame(
            list(range(0, len(psm_report))), 
            columns=['index'], 
            index=psm_report.index
        )), 
        axis=1
    ).fillna('NA')
    
    data = data \
        .groupby(
            ['Peptide', 'Mod', 'Charge', 'GlySite']
        ) \
        .apply(find_best_site) \
        .reset_index(drop=True)
        
    return psm_report.iloc[data['index'].sort_values()]


def extract_assays_from_spectra(psm_report, spectra, 
                                total_fdr_cutoff=0.01,
                                clean_glycan_struct=False,
                                clean_glycan_site=False,
                                return_generator=False):
    assay_builder = GlycoAssayBuilder()
    annotator = SpectrumAnnotator(
        assay_builder=assay_builder
    )
    pglyco = pGlycoToAssayConverter()        
    
    if total_fdr_cutoff is not None:
        psm_report = psm_report \
            .loc[psm_report['TotalFDR'] <= total_fdr_cutoff, :]
    
    if clean_glycan_struct:
        psm_report = remove_suspicious_glycan_struct(psm_report)
        
    if clean_glycan_site:
        psm_report = remove_suspicious_glycan_site(psm_report)
    
    def convert(sp):
        row = np.where(psm_report['PepSpec'] == sp['metadata']['title'])[0]
        if len(row) == 0:
            return None
        
        row = int(row[0])
        info = pglyco.parse_psm_info(psm_report.iloc[row, :])
        sequence = info['peptideSequence']
        glycan_struct = info['glycanStruct']
        glycan_site = info['glycanSite']
        modification = info['modification']
        
        spec = annotator.annotate(
            spectrum=sp, 
            sequence=sequence,
            glycan_struct=glycan_struct,
            glycan_site=glycan_site,
            modification=modification
        )                
        spec.update(info)    
        spec = assay_builder.filter_fragments_by_type(spec)        
        return spec
    
    assays = (
        convert(sp)
        for sp in spectra
        if sp is not None
    )
    assays = (x for x in assays if x is not None)
    if not return_generator:
        assays = list(assays)    
    return assays
        

def extract_assays_from_glabel(psm_report, glabel_report, 
                               matched_ion_has_mz=False,
                               total_fdr_cutoff=0.01,
                               clean_glycan_struct=False,
                               clean_glycan_site=False,
                               return_generator=False):
    pglyco = pGlycoToAssayConverter(matched_ion_has_mz=matched_ion_has_mz)
    assay_builder = GlycoAssayBuilder()
    
    if total_fdr_cutoff is not None:
        psm_report = psm_report \
            .loc[psm_report['TotalFDR'] <= total_fdr_cutoff, :]
    
    if clean_glycan_struct:
        psm_report = remove_suspicious_glycan_struct(psm_report)
    
    if clean_glycan_site:
        psm_report = remove_suspicious_glycan_site(psm_report)
        
    assays = pglyco.report_to_assays(
        psm_report, glabel_report, 
        return_generator=True
    )
    
    def update(assay):
        assay = assay_builder.filter_fragments_by_type(assay) 
        
        if 'fragmentMZ' not in assay['fragments']:
            assay = assay_builder.update_fragment_mz(assay)
        if 'precursorMZ' not in assay:
            assay = assay_builder.update_precursor_mz(assay)
            
        return assay
    
    assays = (update(x) for x in assays)
    
    if not return_generator:
        assays = list(assays)    
    return assays
            
            