from assay.modseq import ModifiedSequenceConverter

def OpenSWATH_columns():        
    seq_converter = ModifiedSequenceConverter()
    def modified_sequence(assay, **kwargs):
        return seq_converter.to_tpp_format(
            sequence=assay['peptideSequence'],
            modification=assay.get('modification', None)
        )
        
    def transition_group_id(assay, **kwargs):
        result = modified_sequence(assay, **kwargs) + '_' + \
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
    
    def convert_transition_type(value):
        if value is None:
            return None
        if isinstance(value, list):
            return [1 if x else 0 for x in value]
        if value:
            return 1
        else:
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
        },
        
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
            'function': modified_sequence,
        },
        {
            'name': 'PrecursorCharge',
            'path': 'precursorCharge',
        },
        {
            'name': 'ProductCharge',
            'path': ['fragments', 'fragmentCharge']
        },
        {
            'name': 'FragmentType',
            'path': ['fragments', 'fragmentType']
        },
        {
            'name': 'FragmentSeriesNumber',
            'path': ['fragments', 'fragmentNumber']
        },
        
        {
            'name': 'TransitionGroupId',
            'function': transition_group_id,
        },
        {
            'name': 'TransitionId',
            'function': transition_id,
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
            'path': ['metadata', 'ce']
        },
        {
            'name': 'LabelType',
            'path': 'labelType',
            'default': 'light'
        }
    ]
    return columns