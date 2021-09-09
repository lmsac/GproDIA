from assay.modseq import ModifiedSequenceConverter 
from pepmass.glycomass import GlycanNode

import re
import itertools

def OpenSWATH_glyco_columns(enable_glycoform_uis=False):       
    seq_converter = ModifiedSequenceConverter()
    
    def modified_sequence(assay, **kwargs):
        return seq_converter.to_tpp_format(
            sequence=assay['peptideSequence'],
            modification=assay.get('modification', None)
        )
    
    def compound_name(assay, **kwargs):
        result = modified_sequence(assay, **kwargs)              
        glycan_struct = assay.get('glycanStruct', None)
        if glycan_struct is not None:
            glycan_site = assay.get('glycanSite', None)  
            if glycan_site is not None:
                result += '_' + str(glycan_site) + ',' + glycan_struct
            else:
                result += '_' + glycan_struct
        metadata = assay.get('metadata', None)
        if metadata is not None:
            if metadata.get('glycanDecoy', False):
                result = 'GLYDECOY_' + result
                decoy = True
            if metadata.get('peptideDecoy', False):
                result = 'PEPDECOY_' + result
                decoy = True
            if metadata.get('decoy', False) and not decoy:
                result = 'DECOY_' + result
                
        return result            
        
    def transition_group_id(assay, **kwargs):
        result = compound_name(assay, **kwargs)        
        result += '_' + \
            str(assay.get('precursorCharge', ''))
        index = kwargs.get('index', None)
        if index is not None:
            result = str(index) + '_' + result
        return result
        
    def transition_id(assay, **kwargs):
        group = transition_group_id(assay, **kwargs)        
        annotation = assay['fragments'].get('fragmentAnnotation', None)
        if annotation is not None:
            result = [
                group + ':' + str(i) + '_' + x
                for i, x in enumerate(annotation)
            ]
        else:
            result = [
                group + ':' + str(i) + '_'
                for i, _ in enumerate(assay['fragments']['fragmentIntensity'])
            ]
        return result   
    
    
    def glycosite_id(assay, **kwargs):
        metadata = assay.get('metadata', None)
        if metadata is not None:
            protein = metadata.get('protein', None)
            protein_site = metadata.get('proteinSite', None)
            if protein is not None and protein_site is not None:
                protein = str(protein).split('/')
                protein_site = str(protein_site).split('/')
                return '/'.join((
                    x + '@' + y 
                    for x, y in zip(protein_site, protein)
                ))
        return None
    
    def glycan_id(assay, **kwargs):
        glycan_struct = assay.get('glycanStruct', None)
        metadata = assay.get('metadata', None)
        if glycan_struct is not None:
            if metadata is not None and metadata.get('glycanDecoy', False):
                result = 'DECOY_' + glycan_struct 
#                precursor_mz = assay.get('precursorMZ', None)
#                if precursor_mz is not None:                    
#                    result += '[' + str(precursor_mz) + ']'
                return result
            else:
                return glycan_struct
        else:
            return None    

    def fragment_series_number(assay, **kwargs):
        fragment_number = assay['fragments'].get('fragmentNumber', None)
        
        glycan_struct = assay.get('glycanStruct', None)
        fragment_glycan = assay['fragments'].get('fragmentGlycan', None)
        if glycan_struct is None or fragment_glycan is None:
            if fragment_number is not None:
                return to_int(fragment_number)
            else:
                return fragment_number
                              
        def parse_series(number, glycan):
            if number is not None:
                return int(number)
            if glycan is None:
                return -1
            if glycan == 'Y0' or glycan == 'Y$':
                return 0
            if re.match('^Y-(([A-Z0-9a-z]+)\\(([0-9]+)\\))+$', glycan) \
                is not None:
                return dict(re.findall('([A-Z0-9a-z]+)\\(([0-9]+)\\)', glycan))                     
            else:
                return -1
        
        if fragment_number is not None:
            series = [
                parse_series(x, fragment_glycan[i])
                for i, x in enumerate(fragment_number)
            ]
        else:
            series = [
                parse_series(None, g)
                for g in fragment_glycan
            ]
            
        monosaccharides = list(set(itertools.chain.from_iterable((
            list(x.keys()) for x in series if isinstance(x, dict)
        ))).union(
            list(GlycanNode.from_str(glycan_struct).composition().keys())
        ))
        monosaccharides.sort()            
            
        return [
            int('0'.join((
                str(x.get(k, 0))
                for k in monosaccharides
            )))
            if isinstance(x, dict) \
            else x
            for x in series
        ]
        

    def to_int(x):
        if x is None:
            return -1
        elif isinstance(x, list):
            return list(map(to_int, x))
        else:
            return int(x)
    
    def convert_transition_type(value):
        if value is None:
            return None
        if isinstance(value, list):
            return [1 if x else 0 for x in value]
        if value:
            return 1
        else:
            return 0
        
    def convert_glycoform(value):
        if value is None:
            return None
        if isinstance(value, list):
            return [
                '|'.join(x) if isinstance(x, list) else x 
                for x in value
            ]
        return value
            
    def to_glycan_composition(x):
        return GlycanNode.from_str(x).composition_str()
    
    if enable_glycoform_uis:
        def transition_id_uis(assay, **kwargs):
            result = transition_id(assay, **kwargs)
            identifying = assay['fragments'].get('identifyingTransition', None)
            glycoform = assay['fragments'].get('glycoform', None)
            decoy = assay['fragments'].get('decoyTransition', None)
            if identifying is not None and glycoform is not None:
                for i, x in enumerate(identifying):
                    if x and glycoform[i] is not None:
                        if decoy is not None and decoy[i]:
                            result[i] += '_UISDECOY_{' 
                        else:
                            result[i] += '_UIS_{' 
                        result[i] += '|'.join(glycoform[i]) + '}'
            return result
        
        def decoy_transition(assay, **kwargs):
            decoy = assay['fragments'].get('decoyTransition', None)
            if decoy is not None:
                return convert_transition_type(decoy)
            
            metadata = assay.get('metadata', None)
            if metadata is not None:
                return 1 if metadata.get('decoy', False) else 0
            return 0
    
    columns = [
        {
            'name': 'PrecursorMz',
            'path': 'precursorMZ'
        },
        {
            'name': 'ProductMz',
            'path': ['fragments', 'fragmentMZ']
        },
        {
            'name': 'LibraryIntensity',
            'path': ['fragments', 'fragmentIntensity']
        },
        {
            'name': 'NormalizedRetentionTime',
            'path': 'rt'
        }
    ]
    

    columns.extend([
        {
            'name': 'GlycoSiteId',
            'function': glycosite_id
        },
        {
            'name': 'ProteinGlycoSite',
            'path': ['metadata', 'proteinSite']
        }        
    ])

    columns.extend([
        {
            'name': 'ProteinId',
            'path': ['metadata', 'protein']
        },
        {
            'name': 'PeptideSequence',
            'path': 'peptideSequence'
        },
        {
            'name': 'ModifiedPeptideSequence',
            'function': modified_sequence
        }
    ])
    
    columns.extend([
        {
            'name': 'GlycanId',
            'function': glycan_id
        },
        {
            'name': 'GlycanStruct',
            'path': 'glycanStruct'
        },
        {
            'name': 'GlycanComposition',
            'path': 'glycanStruct',
            'convert': to_glycan_composition
        },
        {
            'name': 'GlycanSite',
            'path': 'glycanSite'
        },
        {
            'name': 'GlycoPeptideId',
            'function': compound_name
        },
        {
            'name': 'DecoyPeptide',
            'path': ['metadata', 'peptideDecoy'],
            'convert': lambda x: 1 if x else 0,
            'default': 0
        },
        {
            'name': 'DecoyGlycan',
            'path': ['metadata', 'glycanDecoy'],
            'convert': lambda x: 1 if x else 0,
            'default': 0
        },
    ])
     
    columns.extend([
        {
            'name': 'PrecursorCharge',
            'path': 'precursorCharge',
            'convert': to_int
        },
        {
            'name': 'ProductCharge',
            'path': ['fragments', 'fragmentCharge'],
            'convert': to_int
        },
        {
            'name': 'FragmentType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'name': 'FragmentSeriesNumber',
            'function': fragment_series_number
        },
        
        {
            'name': 'TransitionGroupId',
            'function': transition_group_id
        },
        {
            'name': 'TransitionId',
            'function': \
                transition_id_uis \
                if enable_glycoform_uis \
                else transition_id
        },
        {
            'name': 'Decoy',
            'path': ['metadata', 'decoy'],
            'convert': lambda x: 1 if x else 0,
            'default': 0
        },
        {
            'name': 'PeptideGroupLabel',
            'function': transition_group_id
        },
        
        {
            'name': 'DetectingTransition',
            'path': ['fragments', 'detectingTransition'],
            'convert': convert_transition_type,
            'default': 1
        },
        {
            'name': 'IdentifyingTransition',
            'path': ['fragments', 'identifyingTransition'],
            'convert': convert_transition_type,
            'default': 0
        },
        {
            'name': 'QuantifyingTransition',
            'path': ['fragments', 'quantifyingTransition'],
            'convert': convert_transition_type,
            'default': 1
        },        
        
        {
            'name': 'Annotation',
            'path': ['fragments', 'fragmentAnnotation']
        },
        {
            'name': 'CollisionEnergy',
            'path': ['metadata', 'collisionEnergy']
        },
        {
            'name': 'LabelType',
            'path': 'labelType',
            'default': 'light'
        }
    ])
    
    if enable_glycoform_uis:
        columns.extend([
            {
                'name': 'Glycoform',
                'path': ['fragments', 'glycoform'],
                'convert': convert_glycoform
            },
            {
                'name': 'DecoyTransition',
                'path': ['fragments', 'decoyTransition'],
                'function': decoy_transition
            },
        ])
        
    return columns
