from assay.modseq import ModifiedSequenceConverter
from pepmass.glycomass import GlycanNode

def traml_writer_glyco_parameters(): 
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
    
    def fragment_type_name(assay, **kwargs):
        fragment_type = assay['fragments'].get('fragmentType', None)
        fragment_loss = assay['fragments'].get('fragmentLossType', None)
        
        if fragment_type is None:
            return 'non-identified ion'
        
        def to_name(frgtype, frgloss):
            if frgtype not in {'b', 'y', 'a', 'x', 'c', 'z'}:
                return 'non-identified ion'
            if frgloss is None or frgloss == '' or \
                frgloss == 'noloss' or frgloss == 'None':
                return 'frag: ' + frgtype + 'ion'
            elif frgloss in {'NH3', 'H2O'}:
                return 'frag: ' + frgtype + 'ion - ' + frgloss 
            else:
                return 'frag: ' + frgtype + 'ion'
        
        if fragment_loss is not None:
            return [
                to_name(x, fragment_loss[i])
                for i, x in enumerate(fragment_type)
            ]
        else:
            return [
                to_name(x, None)
                for i, x in enumerate(fragment_type)
            ]
      
        
    def convert_decoy_name(decoy):
        if decoy:
            return 'decoy SRM transition'
        else:
            return 'target SRM transition'
        
    def convert_peptide_decoy_name(decoy):
        if decoy:
            return 'peptide_decoy_transition'
        else:
            return None
    
    def convert_glycan_decoy_name(decoy):
        if decoy:
            return 'glycan_decoy_transition'
        else:
            return None
        
            
    def to_glycan_composition(x):
        return GlycanNode.from_str(x).composition_str()
           
    return {
        'protein': {
            'id': {                
                'path': ['metadata', 'protein']
            },
            'protein_accession': {                
                'path': ['metadata', 'protein']
            },
            'sequence': {
                'path': ['metadata', 'proteinSequence']
            }
        },
        'peptide': {
            'id': {                
                'function': transition_group_id
            },
            'sequence': {
                'path': 'peptideSequence'
            },
            'charge_state': {
                'path': 'precursorCharge'
            },            
            'rt': {
                'path': 'rt'
            },                    
            'rt_unit': {
                'path': 'rtUnit'
            },
            'params': [
                {
                    'name': 'peptide group label',
                    'value': {
                        'function': transition_group_id
                    }
                },
                {        
                    'name': 'full_peptide_name',
                    'value': {
                        'function': modified_sequence
                    }
                },
                {
                    'name': 'glycan_struct',
                    'value': {
                        'path': 'glycanStruct'
                    }
                },
                {
                    'name': 'glycan_composition',
                    'value': {
                        'path': 'glycanStruct',
                        'convert': to_glycan_composition
                    }
                },
                {
                    'name': 'glycan_site',
                    'value': {
                        'path': 'glycanSite'
                    }
                },
            ],
        },
        'transition': {
            'id': {                
                'function': transition_id
            },
            'precursor_mz': {
                'path': 'precursorMZ'
            },
            'product_mz': {
                'path': ['fragments', 'fragmentMZ']
            },
            'product_intensity': {
                'path': ['fragments', 'fragmentIntensity']
            },
            'product_charge': {
                'path': ['fragments', 'fragmentCharge']
            },
            
            'interpretation': [
                {
                    'name': 'product ion series ordinal',
                    'value': {
                        'path': ['fragments', 'fragmentNumber']
                    }
                },
                {
                    'name': 'product interpretation rank',
                    'value': 1
                },
                {                        
                    'name': {
                        'function': fragment_type_name
                    }
                },       
            ],         
            
            'params': [                
                {
                    'name': {
                        'path': ['metadata', 'decoy'],
                        'convert': convert_decoy_name
                    }
                },
                {
                    'name': {
                        'path': ['metadata', 'peptideDecoy'],
                        'convert': convert_peptide_decoy_name
                    }
                },
                {
                    'name': {
                        'path': ['metadata', 'glycanDecoy'],
                        'convert': convert_glycan_decoy_name
                    }
                },
                {
                    'name': 'annotation', 
                    'value': {
                        'path': ['fragments', 'fragmentAnnotation']
                    }
                }
            ],
        }
    }