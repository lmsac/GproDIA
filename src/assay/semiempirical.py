from sklearn.neighbors import NearestNeighbors
from collections import Counter
import numpy as np
import pandas as pd
import itertools
import copy

from pepmass.glycomass import GlycanNode
from assay.glycoassay import GlycoAssayBuilder
from assay.consensus import ConsensusAssayCombiner
from assay.assay2table import AssayToDataFrameConverter
from assay.modseq import stringify_modification

def vectorize_glycan(glycan_struct, monosaccharides=None):
    def get_composition(glycan_struct):
        if glycan_struct is None:
            return None
        if isinstance(glycan_struct, str):
            glycan_struct = GlycanNode.from_str(glycan_struct)
        if isinstance(glycan_struct, GlycanNode):
            return glycan_struct.composition()
        else:
            raise TypeError('invalid glycan_struct type:' + type(glycan_struct))
    
    def to_vector(composition, monosaccharides):
        if composition is None:
            if monosaccharides is None:
                return np.zeros(0, dtype=int)
            else:
                return np.zeros(len(monosaccharides), dtype=int)
        else:
            if monosaccharides is None:
                return np.array(list(composition.values()), dtype=int)
            else:
                return np.array(
                    [composition.get(k, 0) for k in monosaccharides], 
                    dtype=int
                )
    
    if glycan_struct is None or \
        isinstance(glycan_struct, str) or \
        isinstance(glycan_struct, GlycanNode):
        composition = get_composition(glycan_struct)
        return to_vector(composition, monosaccharides)
    
    else:
        composition = [
            get_composition(x)
            for x in glycan_struct
        ]        
        if monosaccharides is None:
            monosaccharides = list(set(itertools.chain.from_iterable((
                x.keys()
                for x in composition
                if x is not None
            ))))                    
        return np.array([
            to_vector(x, monosaccharides)
            for x in composition
        ])

def vectorize_peptide(sequence, amino_acids=None):
    def get_composition(sequence):
        if sequence is None:
            return None
        if isinstance(sequence, str):
            return Counter(sequence)
        else:
            raise TypeError('invalid sequence type:' + type(sequence))
    
    def to_vector(composition, amino_acids):
        if composition is None:
            if amino_acids is None:
                return np.zeros(0, dtype=int)
            else:
                return np.zeros(len(amino_acids), dtype=int)
        else:
            if amino_acids is None:
                return np.array(list(composition.values()), dtype=int)
            else:
                return np.array(
                    [composition.get(k, 0) for k in amino_acids], 
                    dtype=int
                )
    
    if sequence is None or \
        isinstance(sequence, str):
        composition = get_composition(sequence)
        return to_vector(composition, amino_acids)
    
    else:
        composition = [
            get_composition(x)
            for x in sequence
        ]        
        if amino_acids is None:
            amino_acids = list(set(itertools.chain.from_iterable((
                x.keys()
                for x in composition
                if x is not None
            ))))                    
        return np.array([
            to_vector(x, amino_acids)
            for x in composition
        ])
            
    
class SemiEmpiricalGlycoAssayBuilder():
    def __init__(self, max_peptide_neighbor_number=3,
                 max_glycan_neighbor_number=3,
                 ignore_modification=False, 
                 assay_builder=None,
                 assay_combiner=None):
        if assay_builder is None:
            assay_builder = GlycoAssayBuilder()
        self.assay_builder = assay_builder
                
        if assay_combiner is None:
            assay_combiner = ConsensusAssayCombiner()
        self.assay_combiner = assay_combiner      
    
        self.ignore_modification = ignore_modification
        self.max_peptide_neighbor_number = max_peptide_neighbor_number
        self.max_glycan_neighbor_number = max_glycan_neighbor_number
        
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
                'convert': str,
                'default': 'None'
            },
            {
                'name': 'precursorCharge',
                'path': 'precursorCharge',
                'convert': str,
                'default': 'None'
            }
        ]
        self.data_converter = AssayToDataFrameConverter(columns=columns)
    
    
    def load_empirical_assays(self, assays=None, 
                              peptide_assays=None, 
                              glycan_assays=None):        
        assays_0 = None
        assay_table = None
        if assays is not None:
            assays_0 = assays
            assay_table = self.data_converter \
                .assays_to_dataframe(assays)
            assay_table['use_peptide'] = True
            assay_table['use_glycan'] = True            

        if peptide_assays is not None:       
            assay_table_1 = self.data_converter \
                .assays_to_dataframe(peptide_assays)
            assay_table_1['use_peptide'] = True
            assay_table_1['use_glycan'] = False
            
            if assay_table is None:
                assays_0 = peptide_assays
                assay_table = assay_table_1
            else:
                assays_0 = assays_0 + peptide_assays
                assay_table = pd.concat(
                    [assay_table, assay_table_1], 
                    ignore_index=True
                )

        if glycan_assays is not None:       
            assay_table_1 = self.data_converter \
                .assays_to_dataframe(glycan_assays)
            assay_table_1['use_peptide'] = False
            assay_table_1['use_glycan'] = True
            if assay_table is None:
                assays_0 = glycan_assays
                assay_table = assay_table_1
            else:
                assays_0 = assays_0 + glycan_assays
                assay_table = pd.concat(
                    [assay_table, assay_table_1], 
                    ignore_index=True
                )
        self.assays = assays_0
        self.assay_table = assay_table
        self.peptide_vectors = {}
        self.glycan_vectors = {}
        self.peptide_fragments = {}
        self.glycan_fragments = {}
    
    
    def find_empirical_assays(self, sequence='any', modification='any', 
                              glycan_struct='any', glycan_site='any', 
                              charge='any', 
                              use_peptide='any', use_glycan='any'):
        assay_table = self.assay_table
        index = assay_table.apply(lambda x: True, axis=1)
        if sequence != 'any':
            index &= assay_table['peptideSequence'] == str(sequence)
        if modification != 'any':
            index &= assay_table['modification'] == \
                str(stringify_modification(modification))
        if glycan_struct != 'any':
            index &= assay_table['glycanStruct'] == str(glycan_struct)
        if glycan_site != 'any':
            index &= assay_table['glycanSite'] == str(glycan_site)
        if charge != 'any':
            index &= assay_table['precursorCharge'] == str(charge)
        
        if use_peptide != 'any':
            index &= assay_table['use_peptide'] == bool(use_peptide)
        if use_glycan != 'any':
            index &= assay_table['use_glycan'] == bool(use_glycan)
        
        return np.where(index)[0]
    
    
    def filter_peptide_empirical_assays_by_glycan_distance(
        self, empirical_assay_indexes, 
        sequence=None, charge=None, modification=None, 
        glycan_struct=None, glycan_site=None):
        monosaccharides=list(self.assay_builder \
                             .mass_calculator.monosaccharide.keys())
        
        def get_glycan_vector(i):
            glycan_vec = self.glycan_vectors.get(i, None)
            if glycan_vec is None:
                glycan_vec = vectorize_glycan(
                    self.assay_table.loc[i, 'glycanStruct'],
                    monosaccharides=monosaccharides
                )
                self.glycan_vectors[i] = glycan_vec
            return glycan_vec
        
        if self.max_glycan_neighbor_number is None or \
            len(empirical_assay_indexes) <= self.max_glycan_neighbor_number:
            n_neighbors = len(empirical_assay_indexes)
        else:
            n_neighbors = self.max_glycan_neighbor_number
        
        vec = np.array([
            get_glycan_vector(i)
            for i in empirical_assay_indexes
        ])            
        x = vectorize_glycan([glycan_struct], monosaccharides=monosaccharides)
        nbrs = NearestNeighbors(
            n_neighbors=n_neighbors, 
            algorithm='ball_tree'
        ).fit(vec)
        distances, indexes = nbrs.kneighbors(x, return_distance=True) 
        
        return empirical_assay_indexes[indexes[0]], distances[0]
    
    
    def filter_glycan_empirical_assays_by_peptide_distance(
        self, empirical_assay_indexes, 
        sequence=None, charge=None, modification=None, 
        glycan_struct=None, glycan_site=None):
        amino_acids = list(self.assay_builder \
                           .mass_calculator.aa_residues.keys())
        
        def get_peptide_vector(i):
            peptide_vec = self.peptide_vectors.get(i, None)
            if peptide_vec is None:
                peptide_vec = vectorize_peptide(
                    self.assay_table.loc[i, 'peptideSequence'],
                    amino_acids=amino_acids
                )
                self.peptide_vectors[i] = peptide_vec
            return peptide_vec
        
        if self.max_peptide_neighbor_number is None or \
            len(empirical_assay_indexes) <= self.max_peptide_neighbor_number:
            n_neighbors = len(empirical_assay_indexes)
        else:
            n_neighbors = self.max_peptide_neighbor_number
        
        vec = np.array([
            get_peptide_vector(i)
            for i in empirical_assay_indexes
        ])            
        x = vectorize_peptide([sequence], amino_acids=amino_acids)
        nbrs = NearestNeighbors(
            n_neighbors=n_neighbors, 
            algorithm='ball_tree'
        ).fit(vec)
        distances, indexes = nbrs.kneighbors(x, return_distance=True) 
        
        return empirical_assay_indexes[indexes[0]], distances[0]
    
    
    def calculate_rt(self, sequence, modification=None, 
                     glycan_struct=None):
        peptide_index = self.find_empirical_assays(
            sequence=sequence, 
            modification=modification \
                if not self.ignore_modification else 'any',
            use_peptide=True
        )        
        peptide_index, peptide_distance = self \
            .filter_peptide_empirical_assays_by_glycan_distance(
                peptide_index, 
                sequence=sequence,
                modification=modification,
                glycan_struct=glycan_struct
            )
        if len(peptide_index) == 0:
            raise ValueError('no peptide empirical assay found: ' + \
                             str(sequence) + ', ' + str(modification))
        
        sigma2 = np.power(peptide_distance, 2).mean()   
        if sigma2 == 0:
            peptide_weight = np.ones_like(peptide_distance)
        else:
            peptide_weight = \
                np.exp(-np.power(peptide_distance, 2) / 2 / sigma2) / \
                np.sqrt(2 * np.pi * sigma2)
                
        rt = float(np.average(
            [
                self.assays[i]['rt']
                for i in peptide_index
            ], 
            weights=peptide_weight
        ))        
        return rt
            
    
    def calculate_fragment_intensity(self, sequence, charge, modification=None, 
                            glycan_struct=None, glycan_site=None):
        peptide_index = self.find_empirical_assays(
            sequence=sequence, 
            modification=modification \
                if not self.ignore_modification else 'any',
            glycan_site=glycan_site \
                if glycan_site is not None else 'any',
            charge=charge,
            use_peptide=True
        )
        peptide_index, peptide_distance = self \
            .filter_peptide_empirical_assays_by_glycan_distance(
                peptide_index, 
                sequence=sequence,
                charge=charge,
                modification=modification,
                glycan_struct=glycan_struct,
                glycan_site=glycan_site
            )
        if len(peptide_index) == 0:
            raise ValueError('no peptide empirical assay found: ' + \
                             str(sequence) + ', ' + str(modification) + ', ' + \
                             str(charge))
        
        glycan_index = self.find_empirical_assays(
            glycan_struct=glycan_struct,
            charge=charge,
            use_glycan=True
        )
        glycan_index, glycan_distance = self \
            .filter_glycan_empirical_assays_by_peptide_distance(
                glycan_index,
                sequence=sequence,
                charge=charge,
                modification=modification,
                glycan_struct=glycan_struct,
                glycan_site=glycan_site
            )
        if len(glycan_index) == 0:
            raise ValueError('no glycan empirical assay found: ' + \
                             str(glycan_struct) + ',' +\
                             str(charge))
        
        def get_peptide_fragments(i):
            fragments = self.peptide_fragments.get(i, None)
            if fragments is None:
                fragments = self.assay_builder.filter_fragments_by_type(
                    self.assays[i],
                    fragment_type=self.assay_builder.fragment_types
                )
                self.peptide_fragments[i] = fragments
                        
            return fragments
        
        def get_glycan_fragments(i):
            fragments = self.glycan_fragments.get(i, None)
            if fragments is None:
                fragments = self.assay_builder.filter_fragments_by_type(
                    self.assays[i],
                    fragment_type=self.assay_builder.glycan_fragment_types
                )                    
                self.glycan_fragments[i] = fragments
                           
            return fragments
                
        def calculate_intensity_ratio(i):
            peptide_fragments = get_peptide_fragments(i)
            glycan_fragments = get_glycan_fragments(i)
            pepsum = np.sum(peptide_fragments['fragments']['fragmentIntensity'])
            glysum = np.sum(glycan_fragments['fragments']['fragmentIntensity'])
            return pepsum / (pepsum + glysum)
        
        
        fragments = {
            'fragmentAnnotation': [],
            'fragmentType': [],
            'fragmentNumber': [],
            'fragmentCharge': [],
            'fragmentLossType': [],
            'fragmentGlycan': []
        }
       
        intensity_ratio = np.mean([
            calculate_intensity_ratio(i)
            for i in np.concatenate((peptide_index, glycan_index))
        ])      
        
        if intensity_ratio > 0:
            peptide_fragments = self.assay_combiner.combine_replicates(
                [get_peptide_fragments(i) for i in peptide_index]
            )
            if peptide_fragments is None:
                return None
            
            for k, v in fragments.items():
                x = peptide_fragments['fragments'].get(k, None)
                if x is not None:
                    v.extend(x)     
            n = max((len(v) for k, v in fragments.items()))
            for k, v in fragments.items():
                if len(v) < n:
                    v.extend([None] * (n - len(v)))
            
        if intensity_ratio < 1:
            glycan_fragments = self.assay_combiner.combine_replicates(
                [get_glycan_fragments(i) for i in glycan_index]
            )
            if glycan_fragments is None:
                return None
            
            for k, v in fragments.items():
                x = glycan_fragments['fragments'].get(k, None)
                if x is not None:
                    v.extend(x)            
            n = max((len(v) for k, v in fragments.items()))
            for k, v in fragments.items():
                if len(v) < n:
                    v.extend([None] * (n - len(v)))
                    
        if intensity_ratio == 1:
            fragments['fragmentIntensity'] = \
                peptide_fragments['fragments']['fragmentIntensity']
        elif intensity_ratio == 0:
            fragments['fragmentIntensity'] = \
                glycan_fragments['fragments']['fragmentIntensity']
        else:  
            pepsum = np.sum(peptide_fragments['fragments']['fragmentIntensity'])
            glysum = np.sum(glycan_fragments['fragments']['fragmentIntensity'])
            if pepsum <= glysum:
                fragments['fragmentIntensity'] = \
                    (np.array(peptide_fragments \
                     ['fragments']['fragmentIntensity']) / pepsum * glysum * \
                     intensity_ratio / (1 - intensity_ratio)).tolist() + \
                    glycan_fragments['fragments']['fragmentIntensity']
            else:
                fragments['fragmentIntensity'] = \
                    peptide_fragments['fragments']['fragmentIntensity'] + \
                    (np.array(glycan_fragments \
                     ['fragments']['fragmentIntensity']) / glysum * \
                     pepsum * (1 - intensity_ratio) / intensity_ratio).tolist()
        
        return fragments
        
            
    def assay(self, sequence, charge, modification=None, 
                    glycan_struct=None, glycan_site=None, **kwargs):
        rt = self.calculate_rt(
            sequence=sequence,
            modification=modification,
            glycan_struct=glycan_struct
        )
        
        fragments = self.calculate_fragment_intensity(
            sequence=sequence,
            charge=charge,
            modification=modification,
            glycan_struct=glycan_struct,
            glycan_site=glycan_site
        )
        
        assay = self.assay_builder.assay(
            sequence=sequence,
            charge=charge,
            modification=modification,
            glycan_struct=glycan_struct,
            glycan_site=glycan_site,
            fragments=fragments,
            rt=rt,
            **kwargs
        )
        
        assay = self.assay_builder.update_precursor_mz(assay)
        assay = self.assay_builder.update_fragment_mz(assay)        
        return assay

    

def get_peptide_table(builder, 
                      min_peptide_occurrence,
                      use_glycan_site=True):
    if 'use_peptide' in builder.assay_table.columns:
        peptides = builder.assay_table[[
            'peptideSequence', 
            'modification',             
            'glycanSite',
            'precursorCharge',
            'use_peptide'
        ]].copy()
        peptides.insert(0, 'index', range(0, len(peptides)))
        peptides = peptides.loc[peptides['use_peptide'] == True, :]
    else:
        peptides = builder.assay_table[[
            'peptideSequence', 
            'modification',             
            'glycanSite',
            'precursorCharge'
        ]].copy()
        peptides.insert(0, 'index', range(0, len(peptides)))
        
    peptides = peptides \
        .groupby([
            'peptideSequence', 
            'modification',             
            'glycanSite',
            'precursorCharge'
        ]) \
        .apply(lambda x: x.iloc[0, :].append( \
               pd.Series(len(x), index=['count']))) \
        .reset_index(drop=True)
    if min_peptide_occurrence is not None:
        peptides = peptides.loc[ \
            peptides['count'] >= min_peptide_occurrence, :]
        
    if not use_glycan_site:
        peptides.sort_values(by='count', ascending=False, inplace=True)
        peptides.drop_duplicates(
            [
                'peptideSequence', 
                'modification',
                'precursorCharge'
            ],
            inplace=True
        )
        peptides.sort_values(by='index', inplace=True)
        
    return peptides


def get_glycan_table(builder, 
                     min_glycan_occurrence,
                     use_glycan_struct=True):
    if 'use_glycan' in builder.assay_table.columns:
        glycans = builder.assay_table[[
            'glycanStruct',
            'precursorCharge',
            'use_glycan'
        ]].copy()
        glycans.insert(0, 'index', range(0, len(glycans)))
        glycans = glycans.loc[glycans['use_glycan'] == True, :]
    else:
        glycans = builder.assay_table[[
            'glycanStruct',
            'precursorCharge'
        ]].copy()
        glycans.insert(0, 'index', range(0, len(glycans)))
    
    glycans = glycans \
        .groupby([
            'glycanStruct',
            'precursorCharge'
        ]) \
        .apply(lambda x: x.iloc[0, :].append( \
               pd.Series(len(x), index=['count']))) \
        .reset_index(drop=True)       
    if min_glycan_occurrence is not None:
        glycans = glycans.loc[ \
            glycans['count'] >= min_glycan_occurrence, :] 
        
    if not use_glycan_struct:
        glycans.sort_values(by='count', ascending=False, inplace=True)
        glycans.insert(
            2, 'glycanComposition', 
            glycans['glycanStruct'] \
                .map(lambda x: GlycanNode.from_str(x).composition_str())
        )
        glycans.drop_duplicates(
            ['glycanComposition', 'precursorCharge'], 
            inplace=True
        )
        glycans.sort_values(by='index', inplace=True)
    
    return glycans


def interchange_peptide_glycan(assays, return_generator=False, 
                               return_glycopeptide_table=False,
                               min_peptide_occurrence=3,
                               min_glycan_occurrence=3,
                               use_glycan_struct=True, 
                               use_glycan_site=True,                               
                               top_n_assays_by_occurrence=None,
                               random_select_n_assays=None,
                               **kwargs):
    builder = SemiEmpiricalGlycoAssayBuilder(**kwargs)
    builder.load_empirical_assays(assays)
    
    glycopeptides_existed = builder.assay_table[[
        'peptideSequence', 'modification', 
        'glycanStruct', 'glycanSite', 
        'precursorCharge'
    ]].copy()
    
    if not use_glycan_site:
        glycopeptides_existed.drop(columns=['glycanSite'],inplace=True)
    
    if not use_glycan_struct:        
        glycopeptides_existed.insert(
            3, 'glycanComposition', 
            glycopeptides_existed['glycanStruct'] \
                .map(lambda x: GlycanNode.from_str(x).composition_str())
        )
        glycopeptides_existed.drop(columns=['glycanStruct'],inplace=True)

    glycopeptides_existed.drop_duplicates(inplace=True)
    
    peptides = get_peptide_table(
        builder, 
        min_peptide_occurrence=min_peptide_occurrence,
        use_glycan_site=use_glycan_site
    )
    glycans = get_glycan_table(
        builder,
        min_glycan_occurrence=min_glycan_occurrence,
        use_glycan_struct=use_glycan_struct
    )    
    glycopeptides = pd.merge(
        peptides, glycans, 
        on='precursorCharge', 
        suffixes=['_peptide', '_glycan']
    )
    glycopeptides = glycopeptides \
        .merge(
            glycopeptides_existed, 
            how='outer', indicator=True
        ) \
        .loc[lambda x : x['_merge']=='left_only'] \
        .drop(columns='_merge')
    
    if top_n_assays_by_occurrence is not None and \
        len(glycopeptides) > top_n_assays_by_occurrence:
        glycopeptides = glycopeptides.iloc[np.argsort( \
            -(glycopeptides['count_peptide'] + glycopeptides['count_glycan']) \
        )[:top_n_assays_by_occurrence], :]
            
    if random_select_n_assays is not None and \
        len(glycopeptides) > random_select_n_assays:
        glycopeptides = glycopeptides.sample(random_select_n_assays)
        
    assays = (
        builder.assay(
            sequence=builder.assays[i]['peptideSequence'], 
            charge=builder.assays[i]['precursorCharge'], 
            modification=builder.assays[i].get('modification', None), 
            glycan_struct=builder.assays[j]['glycanStruct'], 
            glycan_site=builder.assays[i].get('glycanSite', None), 
            metadata=copy.deepcopy(builder.assays[i].get('metadata', None))
        )
        for i, j in zip(
            glycopeptides['index_peptide'].astype(int), 
            glycopeptides['index_glycan'].astype(int)
        )
    )
    assays = (x for x in assays if x is not None)
    if not return_generator:
        assays = list(assays)
        
    if return_glycopeptide_table:
        return assays, glycopeptides
    else:
        return assays


def exchange_peptide_glycan(peptide_assays, glycan_assays, 
                            return_generator=False,
                            return_glycopeptide_table=False,
                            min_peptide_occurrence=3,
                            min_glycan_occurrence=3,
                            use_glycan_struct=True,
                            use_glycan_site=True,
                            top_n_assays_by_occurrence=None,
                            random_select_n_assays=None,
                            **kwargs):
    builder = SemiEmpiricalGlycoAssayBuilder(**kwargs)
    builder.load_empirical_assays(
        peptide_assays=peptide_assays,
        glycan_assays=glycan_assays
    )
    
    glycopeptides_existed = builder.assay_table[[
        'peptideSequence', 'modification', 
        'glycanStruct', 'glycanSite', 
        'precursorCharge'
    ]].copy()
    
    if not use_glycan_site:
        glycopeptides_existed.drop(columns=['glycanSite'],inplace=True)
    
    if not use_glycan_struct:        
        glycopeptides_existed.insert(
            3, 'glycanComposition', 
            glycopeptides_existed['glycanStruct'] \
                .map(lambda x: GlycanNode.from_str(x).composition_str())
        )
        glycopeptides_existed.drop(columns=['glycanStruct'],inplace=True)

    glycopeptides_existed.drop_duplicates(inplace=True)
    
    peptides = get_peptide_table(
        builder, 
        min_peptide_occurrence=min_peptide_occurrence,
        use_glycan_site=use_glycan_site
    )
    glycans = get_glycan_table(
        builder,
        min_glycan_occurrence=min_glycan_occurrence,
        use_glycan_struct=use_glycan_struct
    )
    
    glycopeptides = pd.merge(
        peptides, glycans, 
        on='precursorCharge', 
        suffixes=['_peptide', '_glycan']
    )
    glycopeptides = glycopeptides \
        .drop(columns=['use_peptide', 'use_glycan']) \
        .merge(
            glycopeptides_existed, 
            how='outer', indicator=True
        ) \
        .loc[lambda x : x['_merge']=='left_only'] \
        .drop(columns='_merge')
        
    if top_n_assays_by_occurrence is not None and \
        len(glycopeptides) > top_n_assays_by_occurrence:
        glycopeptides = glycopeptides.iloc[np.argsort( \
            -(glycopeptides['count_peptide'] + glycopeptides['count_glycan']) \
        )[:top_n_assays_by_occurrence], :]
    
    if random_select_n_assays is not None and \
        len(glycopeptides) > random_select_n_assays:
        glycopeptides = glycopeptides.sample(random_select_n_assays)
    
    assays = (
        builder.assay(
            sequence=builder.assays[i]['peptideSequence'], 
            charge=builder.assays[i]['precursorCharge'], 
            modification=builder.assays[i].get('modification', None), 
            glycan_struct=builder.assays[j]['glycanStruct'], 
            glycan_site=builder.assays[i].get('glycanSite', None), 
            metadata=copy.deepcopy(builder.assays[i].get('metadata', None))
        )
        for i, j in zip(
            glycopeptides['index_peptide'].astype(int), 
            glycopeptides['index_glycan'].astype(int)
        )
    )
    assays = (x for x in assays if x is not None)
    if not return_generator:
        assays = list(assays)
        
    if return_glycopeptide_table:
        return assays, glycopeptides
    else:
        return assays


def generate_semiempirical_assays_for_cross_validation(
    assays, include_index=False, 
    return_generator=False, 
    return_glycopeptide_table=False, 
    min_peptide_occurrence=3, min_glycan_occurrence=3,
    random_select_n_assays=None,
    **kwargs):
    
    builder = SemiEmpiricalGlycoAssayBuilder(**kwargs)    
    builder.load_empirical_assays(assays)
    
    glycopeptides = builder.assay_table[[
        'peptideSequence', 
        'modification',             
        'glycanSite',
        'glycanStruct',
        'precursorCharge'
    ]].copy()
    glycopeptides.insert(0, 'index', range(0, len(glycopeptides)))
        
    glycopeptides = glycopeptides \
        .groupby([
            'peptideSequence', 
            'modification',             
            'glycanSite',
            'precursorCharge'
        ]) \
        .apply(lambda x: pd.concat(
            (x, pd.DataFrame(
                [len(x)] * len(x), 
                columns=['count_peptide'], 
                index=x.index
            )),
            axis=1
        )) \
        .reset_index(drop=True)
        
    glycopeptides = glycopeptides \
        .groupby([
            'glycanStruct',
            'precursorCharge'
        ]) \
        .apply(lambda x: pd.concat(
            (x, pd.DataFrame(
                [len(x)] * len(x), 
                columns=['count_glycan'], 
                index=x.index
            )),
            axis=1
        )) \
        .reset_index(drop=True)
    
    if min_peptide_occurrence is None:
        min_peptide_occurrence = 1
    glycopeptides = glycopeptides.loc[ \
        glycopeptides['count_peptide'] >= min_peptide_occurrence + 1, :]  
    
    if min_glycan_occurrence is None:
        min_glycan_occurrence = 1
    glycopeptides = glycopeptides.loc[ \
        glycopeptides['count_glycan'] >= min_glycan_occurrence + 1, :]    
    
    if random_select_n_assays is not None and \
        len(glycopeptides) > random_select_n_assays:
        glycopeptides = glycopeptides.sample(random_select_n_assays)
    
    
    assay_table = builder.assay_table
    peptide_vectors = {}
    glycan_vectors = {}
    peptide_fragments = {}
    glycan_fragments = {}
    
    def generate_assay(i):
        builder.assay_table = assay_table.drop(i, axis=0) \
            .reset_index(drop=True)
        builder.assays = [x for j, x in enumerate(assays) if j != i]    
        
        def create_builder_cache(cache):
            return {
                (j if j < i else j - 1): v 
                for j, v in cache.items() 
                if j != i
            }
        
        builder.peptide_vectors = create_builder_cache(peptide_vectors)
        builder.glycan_vectors = create_builder_cache(glycan_vectors)
        builder.peptide_fragments = create_builder_cache(peptide_fragments)
        builder.glycan_fragments = create_builder_cache(glycan_fragments)
        
        # builder = SemiEmpiricalGlycoAssayBuilder(**kwargs)
        # assays_other = [x for j, x in enumerate(assays) if j != i]
        # builder.load_empirical_assays(assays_other)
        
        result = builder.assay(
            sequence=assays[i]['peptideSequence'], 
            charge=assays[i]['precursorCharge'], 
            modification=assays[i].get('modification', None), 
            glycan_struct=assays[i]['glycanStruct'], 
            glycan_site=assays[i].get('glycanSite', None), 
            metadata=copy.deepcopy(assays[i].get('metadata', None))
        )
        
        def update_cache(cache, builder_cache):
            for j, v in builder_cache.items():
                if j >= i:
                    j += 1
                cache.setdefault(j, v)
        
        update_cache(peptide_vectors, builder.peptide_vectors)
        update_cache(glycan_vectors, builder.glycan_vectors)
        update_cache(peptide_fragments, builder.peptide_fragments)
        update_cache(glycan_fragments, builder.glycan_fragments)
        
        return result
        
    new_assays = (
        (i, generate_assay(i))
        for i in glycopeptides['index'].astype(int)
    )
    new_assays = (
        x if include_index else x[1] 
        for x in new_assays 
        if x[1] is not None
    )
    if not return_generator:
        new_assays = list(new_assays)
        
    if return_glycopeptide_table:
        return new_assays, glycopeptides
    else:
        return new_assays


def generate_semiempirical_assays(data, assays, 
                                  peptide_assays=None, glycan_assays=None,
                                  return_generator=False,
                                  return_glycopeptide_table=False,
                                  min_peptide_occurrence=3,
                                  min_glycan_occurrence=3,
                                  use_glycan_struct=True,
                                  use_glycan_site=True,
                                  **kwargs):
    builder = SemiEmpiricalGlycoAssayBuilder(**kwargs)
    builder.load_empirical_assays(
        assays=assays,
        peptide_assays=peptide_assays,
        glycan_assays=glycan_assays
    )
    
    glycopeptides = data.copy()
    if not use_glycan_site and 'glycanSite' in glycopeptides.columns:
        glycopeptides.drop(columns=['glycanSite'], inplace=True)
    if not use_glycan_struct and 'glycanStruct' in glycopeptides.columns:
        if not 'glycanComposition' in glycopeptides.columns:
            glycopeptides['glycanComposition'] = \
                glycopeptides['glycanStruct'] \
                    .map(lambda x: GlycanNode.from_str(x).composition_str())
        glycopeptides.drop(columns=['glycanStruct'], inplace=True)
    
    peptides = get_peptide_table(
        builder, 
        min_peptide_occurrence=min_peptide_occurrence,
        use_glycan_site=use_glycan_site
    )
    glycans = get_glycan_table(
        builder,
        min_glycan_occurrence=min_glycan_occurrence,
        use_glycan_struct=use_glycan_struct
    )

    glycopeptides = pd.merge(
        glycopeptides, 
        peptides.rename(columns={
            k: k + '_peptide'
            for k in ['index', 'count']
        })
    )    
    glycopeptides = pd.merge(
        glycopeptides, 
        glycans.rename(columns={
            k: k + '_glycan'
            for k in ['index', 'count']
        })
    )
    
    if 'use_peptide' in glycopeptides.columns:
        glycopeptides.drop(columns=['use_peptide'], inplace=True)
    if 'use_glycan' in glycopeptides.columns:
        glycopeptides.drop(columns=['use_glycan'], inplace=True)
    
    assays = (
        builder.assay(
            sequence=builder.assays[i]['peptideSequence'], 
            charge=builder.assays[i]['precursorCharge'], 
            modification=builder.assays[i].get('modification', None), 
            glycan_struct=builder.assays[j]['glycanStruct'], 
            glycan_site=builder.assays[i].get('glycanSite', None), 
            metadata=copy.deepcopy(builder.assays[i].get('metadata', None))
        )
        for i, j in zip(
            glycopeptides['index_peptide'].astype(int), 
            glycopeptides['index_glycan'].astype(int)
        )
    )
    assays = (x for x in assays if x is not None)
    if not return_generator:
        assays = list(assays)
        
    if return_glycopeptide_table:
        return assays, glycopeptides
    else:
        return assays
    
    