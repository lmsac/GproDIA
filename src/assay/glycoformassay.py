import copy
import itertools
import numpy as np
import pandas as pd
import random

from .glycoassay import GlycoAssayBuilder
from .assay2table import AssayToDataFrameConverter
from .modseq import stringify_modification


class GlycoformUisAssayBuilder(GlycoAssayBuilder):
    def __init__(self,
                 background_glycans=None,
                 max_background_glycan_number=10,
                 enable_identification_ms2_precursors=True,
                 **kwargs):
        super(GlycoformUisAssayBuilder, self).__init__(
            **kwargs
        )

        columns = [
            {
                'name': 'peptideSequence',
                'path': 'peptideSequence',
                'convert': str,
                'default': 'None'
            },
            {
                'name': 'modification',
                'path': 'modification',
                'convert': lambda x: str(stringify_modification(x)),
                'default': 'None'
            },
            {
                'name': 'glycanStruct',
                'path': 'glycanStruct',
                'convert': str,
                'default': 'None'
            },
            {
                'name': 'glycanSite',
                'path': 'glycanSite',
                'default': -1
            },
            {
                'name': 'precursorCharge',
                'path': 'precursorCharge',
                'default': 0
            },
            {
                'name': 'precursorMZ',
                'path': 'precursorMZ',
                'default': -1
            }
        ]
        self.data_converter = AssayToDataFrameConverter(columns=columns)

        if background_glycans is None:
            self.background_glycans = None
        else:
            self.background_glycans = \
                pd.DataFrame(background_glycans, columns=['glycanStruct'])
            self.background_glycans['mw'] = \
                self.background_glycans['glycanStruct'] \
                    .map(lambda x: self.mass_calculator.glycan_mw(x))

        self.max_background_glycan_number = max_background_glycan_number
        self.enable_identification_ms2_precursors = \
            enable_identification_ms2_precursors


    def group_assays(self, assays, swath_windows=None,
                     return_index=False):
        data = self.data_converter.assays_to_dataframe(assays)
        data['index'] = np.arange(0, len(data), dtype=int)

        if swath_windows is not None:
            data['isolationWindow'] = data['precursorMZ'].map(lambda x: \
                np.max(np.append(np.where(
                    (swath_windows['start'] < x) & (swath_windows['end'] > x)
                )[0], -1))
            )
            indexes = list(data \
                .groupby(by=[
                    'peptideSequence',
                    'modification',
                    'isolationWindow'
                ]) \
                .apply(lambda x: x['index'].values.tolist())
            )
        else:
            indexes = list(data \
                .groupby(by=[
                    'peptideSequence',
                    'modification'
                ]) \
                .apply(lambda x: x['index'].values.tolist())
            )
        if return_index:
            return indexes
        else:
            return ([assays[i] for i in x] for x in indexes)


    def get_background_glycans(self, assay, swath_windows=None):
        if self.background_glycans is None:
            return None

        precursor_mz = assay['precursorMZ']

        if swath_windows is not None:
            isolation_window = np.max(np.append(np.where(
                (swath_windows['start'] < precursor_mz) & \
                (swath_windows['end'] > precursor_mz)
            )[0], -1))
            if isolation_window < 0:
                return None
            glycan_mw = self.mass_calculator.glycan_mw(assay['glycanStruct'])

            background_mz = precursor_mz + \
                (self.background_glycans['mw'] - glycan_mw) / \
                assay['precursorCharge']

            return list(self.background_glycans.loc[ \
                (background_mz > swath_windows['start'][isolation_window]) & \
                (background_mz < swath_windows['end'][isolation_window]) \
            ]['glycanStruct'])

        else:
            return list(self.background_glycans['glycanStruct'])


    def filter_background_glycans(self, background_glycans, fragment_names):
        if background_glycans is None or len(background_glycans) == 0:
            return background_glycans

        if not hasattr(self, 'background_glycan_fragments'):
            self.background_glycan_fragments = {}

        all_fragment_names = set(itertools.chain.from_iterable(fragment_names))
        fragment_names = fragment_names.copy()
        exclude = []
        jaccard = []
        for glycan in background_glycans:
            fragments = self.background_glycan_fragments.get(glycan, None)
            if fragments is None:
                fragments = [
                    'Y-' + ''.join(
                        k + '(' + str(v) + ')'
                        for k, v in x.items()
                    )
                    for x in self.mass_calculator.glycan_fragment(glycan)
                ]
                self.background_glycan_fragments[glycan] = fragments

            if any((
                set(frag_names) == set(fragments)
                for frag_names in fragment_names
            )):
                exclude.append(True)
                continue

            exclude.append(False)
            fragment_names.append(fragments)
            jaccard.append(
                len(set(fragments).intersection(all_fragment_names)) / \
                len(set(fragments).union(all_fragment_names))
            )

        background_glycans = [
            x for i, x in enumerate(background_glycans)
            if not exclude[i]
        ]
        indexes = np.argsort(jaccard)[::-1]
        if self.max_background_glycan_number is not None:
            indexes = indexes[:self.max_background_glycan_number]
        return [background_glycans[i] for i in indexes]


    def generate_identifying_transitions(self, assay_group,
                                         swath_windows=None,
                                         background_glycans=None):
        def theoretical_assay(sequence, modification,
                              glycan_struct, glycan_site,
                              precursor_charge):
            assay = self.theoretical_fragments(
                sequence=sequence,
                modification=modification,
                glycan_struct=glycan_struct,
                glycan_site=glycan_site,
                fragment_type=self.glycan_fragment_types
            )
            assay = self.filter_glycan_fragments_by_monosaccharide_number(
                assay,
                min_monosaccharide_number=1
            )   
            assay['precursorCharge'] = precursor_charge
            assay = self.update_precursor_mz(assay)
            if swath_windows is not None:
                assay = self.exclude_fragments_in_isolation_window(
                    assay,
                    swath_windows=swath_windows
                )
            return assay

        theoretical_assays = [
            theoretical_assay(
                sequence=assay['peptideSequence'],
                modification=assay.get('modification', None),
                glycan_struct=assay.get('glycanStruct', None),
                glycan_site=assay.get('glycanSite', None),
                precursor_charge=assay['precursorCharge']
            )
            for assay in assay_group
        ]

        frag_dict = dict()
        for i, assay_1 in enumerate(theoretical_assays):
            for x in assay_1['fragments']['fragmentGlycan']:
                l = frag_dict.get(x, None)
                if l is None:
                    l = set()
                    frag_dict[x] = l
                l.add(i)

        if background_glycans is not None and len(background_glycans) > 0:
            background_glycans = list(set(background_glycans).difference([
                x.get('glycanStruct', None)
                for x in theoretical_assays
            ]))
            if len(background_glycans) > 0:
                background_glycans = self.filter_background_glycans(
                    background_glycans,
                    fragment_names=[
                        x['fragments']['fragmentGlycan']
                        for x in theoretical_assays
                    ]
                )

        if background_glycans is not None and len(background_glycans) > 0:
            background_assays = [
                theoretical_assay(
                    sequence=assay_group[0]['peptideSequence'],
                    modification=assay_group[0].get('modification', None),
                    glycan_struct=glycan,
                    glycan_site=assay_group[0].get('glycanSite', None),
                    precursor_charge=assay_group[0]['precursorCharge']
                )
                for glycan in background_glycans
            ]
            for i, assay_1 in enumerate(background_assays):
                for x in assay_1['fragments']['fragmentGlycan']:
                    l = frag_dict.get(x, None)
                    if l is None:
                        l = set()
                        frag_dict[x] = l
                    l.add(i + len(theoretical_assays))
            theoretical_assays.extend(background_assays)            
   
        identifying_transitions = {
            'glycoform': [],
            'fragmentGlycan': [],
            'fragmentMZ': [],
            'fragmentType':[],
            'fragmentCharge': [],
            'fragmentAnnotation': []
        }
        
        for k, v in frag_dict.items():
            v = list(v)
            for i, x in enumerate(theoretical_assays[v[0]] \
                                  ['fragments']['fragmentGlycan']):
                if k != x:
                    continue
                
                for k1, v1 in theoretical_assays[v[0]]['fragments'].items():
                    l = identifying_transitions.get(k1, None)
                    if l is None:
                        l = [None] * \
                            len(identifying_transitions['glycoform'])
                        identifying_transitions[k1] = l
                    l.append(theoretical_assays[v[0]]['fragments'][k1][i])

                identifying_transitions['glycoform'].append([
                    theoretical_assays[j]['glycanStruct']
                    for j in v
                ])                

        if self.enable_identification_ms2_precursors:
            prec_dict = {}
            for i, assay_1 in enumerate(theoretical_assays):
                prec_mz = assay_1['precursorMZ']
                l = prec_dict.get(prec_mz, None)
                if l is None:
                    l = set()
                    prec_dict[prec_mz] = l
                l.add(i)

            for k, v in prec_dict.items():
                v = list(v)
                for k1, v1 in identifying_transitions.items():
                    if k1 == 'fragmentMZ':
                        v1.append(k)
                    elif k1 == 'fragmentType':
                        v1.append('')
                    elif k1 == 'fragmentCharge':
                        v1.append(theoretical_assays[v[0]]['precursorCharge'])
                    elif k1 == 'fragmentAnnotation':
                        v1.append('MS2_Precursor_i0')
                    elif k1 == 'glycoform':
                        v1.append([
                            theoretical_assays[j]['glycanStruct']
                            for j in v
                        ])
                    else:
                        v1.append(None)

        identifying_transitions['identifyingTransition'] = \
            [True] * len(identifying_transitions['fragmentGlycan'])
        identifying_transitions['detectingTransition'] = \
            [False] * len(identifying_transitions['fragmentGlycan'])
        identifying_transitions['quantifyingTransition'] = \
            [False] * len(identifying_transitions['fragmentGlycan'])

        result = []
        for i, assay in enumerate(assay_group):
            new_assay = copy.deepcopy(assay)
            new_assay['fragments']['detectingTransition'] = \
                [True] * len(assay['fragments']['fragmentMZ'])
            new_assay['fragments']['identifyingTransition'] = \
                [False] * len(assay['fragments']['fragmentMZ'])

            for k in set(new_assay['fragments'].keys()) \
                .union(identifying_transitions.keys()):
                list1 = new_assay['fragments'].get(k, None)
                list2 = identifying_transitions.get(k, None)
                if not isinstance(list1, list):
                    list1 = [list1] * len(assay['fragments']['fragmentMZ'])
                if not isinstance(list2, list):
                    list2 = [list2] * len(identifying_transitions['fragmentMZ'])
                new_assay['fragments'][k] = list1 + list2

            result.append(new_assay)

        return result


    def generate_decoy_identifying_transitions(self, assay_group):
        result = []

        sequence = assay_group[0]['peptideSequence']
        aa_residues = list(self.mass_calculator.aa_residues.keys())

        for i in range(100):
            decoy_sequence = ''.join((
                random.choice(aa_residues)
                for i in range(len(sequence))
            ))
            decoy_mw_shift = self.mass_calculator.mw(sequence=decoy_sequence) - \
                self.mass_calculator.mw(sequence=sequence)
            if abs(decoy_mw_shift) > 1:
                break
            elif i == 99:
                raise ValueError('decoy sequence max attempts')

        for i, assay in enumerate(assay_group):
            new_assay = copy.deepcopy(assay)
            if new_assay['fragments'].get('decoyTransition', None) is None:
                new_assay['fragments']['decoyTransition'] = \
                    [False] * len(assay['fragments']['fragmentMZ'])


            for i, x in enumerate(assay['fragments']['fragmentMZ']):
                if not new_assay['fragments']['identifyingTransition'][i] or \
                    new_assay['fragments']['decoyTransition'][i]:
                    continue

                decoy_mz_shift = decoy_mw_shift / \
                    new_assay['fragments']['fragmentCharge'][i]
                new_assay['fragments']['decoyTransition'].append(True)
                new_assay['fragments']['fragmentMZ'].append(
                    new_assay['fragments']['fragmentMZ'][i] + \
                    decoy_mz_shift
                )
                new_assay['fragments']['fragmentAnnotation'].append(
                    'DECOY_' + \
                    new_assay['fragments']['fragmentAnnotation'][i] + \
                    '[' + ('+' if decoy_mz_shift > 0 else '') + \
                    str(decoy_mz_shift) + ']'
                )
                for k in new_assay['fragments'].keys():
                    if k not in {
                        'fragmentMZ', 'fragmentAnnotation',
                        'decoyTransition'
                    }:
                        new_assay['fragments'][k].append(
                            new_assay['fragments'][k][i]
                        )

            result.append(new_assay)

        return result


    def build_uis_assays(self, assays, swath_windows=None,
                         min_fragment_mz=None,
                         max_fragment_mz=None,
                         return_generator=False):
        def add_identifying_transitions(assay_group):
            background_glycans = self.get_background_glycans(
                assay_group[0],
                swath_windows=swath_windows
            )

            assay_group = self.generate_identifying_transitions(
                assay_group,
                background_glycans=background_glycans,
                swath_windows=swath_windows,
            )
            assay_group = self.generate_decoy_identifying_transitions(
                assay_group
            )
            
            if min_fragment_mz is not None or max_fragment_mz is not None:
                assay_group = [
                    self.filter_fragments_by_mz(
                        assay, 
                        min_mz=min_fragment_mz,
                        max_mz=max_fragment_mz
                    )
                    for assay in assay_group
                ]
            return assay_group

        result = itertools.chain.from_iterable((
            add_identifying_transitions(assay_group)
            for assay_group in self.group_assays(
                assays,
                swath_windows=swath_windows
            )
        ))
        if not return_generator:
            result = list(result)
        return result

