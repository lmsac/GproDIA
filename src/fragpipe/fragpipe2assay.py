import re
import numpy as np

import pandas as pd

from assay.annotation import SpectrumAnnotator
from assay.glycoassay import GlycoAssayBuilder
from assay.modseq import ModifiedSequenceConverter
from pepmass.glycomass import GlycoPeptideMassCalculator, GlycanNode

class FragPipeToAssayConverter:
    def __init__(self,
                 glycan_struct,
                 modseq_converter=None):
        if modseq_converter is not None:
            self.modseq_converter = modseq_converter
        else:
            self.modseq_converter = ModifiedSequenceConverter(GlycoPeptideMassCalculator())

        self.monosaccharides = {
            'Hex': 'H',
            'HexNAc': 'N',
            'NeuAc': 'A',
            'NeuGc': 'G',
            'Fuc': 'F',
        }

        self.glycan_struct = pd.DataFrame.from_dict({
            'glycan_composition': [
                GlycanNode.from_str(gly).composition_str()
                for gly in glycan_struct
            ],
            'glycan_struct': glycan_struct,
        })


    def parse_modified_sequence(self, modified_sequence):
        modified_sequence = re.sub("N\\[[0-9]+\\]", "J", modified_sequence)
        sequence, modification = self.modseq_converter.parse_tpp_format(modified_sequence)
        glycan_site = sequence.index('J') + 1
        return sequence, modification, glycan_site

    def parse_glycan(self, glycan_composition):
        glycan_composition = glycan_composition.split(" % ")[0]

        try:
            glycan_composition = ''.join(
                m + '(' + str(v) + ')' for m, v in sorted((
                    (self.monosaccharides[m_], v_)
                    for m_, v_ in re.findall("([A-Za-z0-9]+)\\(([0-9]+)\\)", glycan_composition)
                ), key=lambda t: t[0])
            )
        except KeyError as e:
            import warnings
            warnings.warn('monosaccharides not found: ' + str(e.args))
            return None

        index = np.where(self.glycan_struct['glycan_composition'] == glycan_composition)[0]
        if len(index) == 0:
            import warnings
            warnings.warn('glycan not found: ' + glycan_composition)
            return None

        index = int(index[0])
        glycan_struct = self.glycan_struct['glycan_struct'].iloc[index]

        return glycan_struct


    def parse_psm_info(self, psm):
        sequence, modification, glycan_site = self.parse_modified_sequence(psm['Modified Peptide'])
        charge = int(psm['Charge'])

        glycan_struct = self.parse_glycan(psm['Observed Modifications'])

        result = {
            'peptideSequence': sequence,
            'precursorCharge': charge,
            'modification': modification,
            'glycanStruct': glycan_struct,
            'glycanSite': glycan_site
        }

        rt = psm.get('Retention', None)
        if rt is not None:
            result.update({
                'rt': float(rt),
            })

        precursor_mz = psm.get('Observed M/Z', None)
        if precursor_mz is not None:
            result.update({
                'precursorMZ': float(precursor_mz),
            })

        glycan_fdr = psm.get('Glycan q-value', None)
        if glycan_fdr is not None:
            result['metadata'] = {
                'glycanFDR': float(glycan_fdr),
            }

        glycan_score = psm.get('Glycan Score', None)
        peptide_score = psm.get('Hyperscore', None)
        total_score = psm.get('Best Score with Delta Mass', None)
        if total_score is not None and \
            peptide_score is not None and \
            glycan_score is not None:
            result['metadata'].update({
                'glycanScore': float(glycan_score),
                'peptideScore': float(peptide_score),
                'totalScore': float(total_score)
            })

        protein = psm.get('Protein', None)
        protein_site = psm.get('Protein Start', None)
        if protein_site is not None:
            protein_site += glycan_site - 1
        if protein is not None:
            result['metadata'].update({
                'protein': protein,
                'proteinSite': protein_site
            })

        spectrum = psm.get('Spectrum', None)
        if spectrum is not None:
            file = re.sub('\\.[0-9]+\\.[0-9]+\\.[0-9]+$', '', spectrum)
            scan = spectrum.split('.')[-2]
        if file is not None and scan is not None:
            result['metadata'].update({
                'file': file,
                'scan': int(scan)
            })

        return result


def extract_assays_from_spectra(psm_report, spectra,
                                glycan_struct,
                                glycan_fdr_cutoff=0.01,
                                return_generator=False):
    assay_builder = GlycoAssayBuilder()
    annotator = SpectrumAnnotator(
        assay_builder=assay_builder
    )
    fragpipe = FragPipeToAssayConverter(glycan_struct)

    if glycan_fdr_cutoff is not None:
        psm_report = psm_report \
            .loc[psm_report['Glycan q-value'] <= glycan_fdr_cutoff, :]

    psm_report.insert(0, 'SpectrumTitle', psm_report['Spectrum'].str.replace("\\.0*([0-9]+)\\.[0-9]+\\.[0-9]+$", ".\\1.\\1."))

    def convert(sp):
        row = np.where(psm_report['SpectrumTitle'] == sp['metadata'].get('title', '').split(' File:')[0])[0]
        if len(row) == 0:
            return None

        row = int(row[0])
        info = fragpipe.parse_psm_info(psm_report.iloc[row, :])

        sequence = info['peptideSequence']
        glycan_struct = info['glycanStruct']
        glycan_site = info['glycanSite']
        modification = info['modification']

        if glycan_struct is None:
            return None

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





