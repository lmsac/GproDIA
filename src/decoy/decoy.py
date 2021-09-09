import copy

from assay import AssayBuilder

class DecoyAssayGenerator():
    def __init__(self, assay_builder=None, 
                 unknown_fragment_type='keep',
                 **kwargs):        
        if assay_builder is None:
            assay_builder = AssayBuilder()
        self.assay_builder = assay_builder
        
        if unknown_fragment_type not in {'ignore', 'keep'}:
            raise ValueError('unknown_fragment_type should be "ignore", "keep" or "error"')
        self.unknown_fragment_type = unknown_fragment_type
        
    
    def decoy_sequence(self, assay):
        sequence = assay['peptideSequence']
        modification = assay.get('modification', None)
        
        if sequence[-1] in {'R', 'K'}:
            index = list(range(len(sequence) - 2, -1, -1)) + \
                [len(sequence) - 1]
        else:
            index = list(range(len(sequence) - 1, -1, -1))
        
        new_sequence = ''.join(sequence[i] for i in index)
        
        if modification is None:
            new_modification = None
        else:
            new_modification = []            
            for i, x in enumerate(index):
                for m in modification:                    
                    if m['position'] == x + 1:
                        m = m.copy()
                        m['position'] = i + 1
                        new_modification.append(m)
            for m in modification:
                if m['position'] is None:                    
                    if m['site'] == 'N-term':
                        new_modification.insert(0, m.copy())
                    elif m['site'] == 'C-term':
                        new_modification.append(m.copy())
        
        return {
            'peptideSequence': new_sequence,
            'modification': new_modification
        }
    
        
    def update_fragments(self, assay):
        fragments = copy.deepcopy(assay['fragments'])
        
        if self.unknown_fragment_type != 'error':
            fragment_index = self.assay_builder.filter_fragments_by_type(
                assay, return_index=True
            )
            assay = self.assay_builder.filter_fragments_by_index(
                assay, fragment_index=fragment_index
            )
        
        assay = self.assay_builder.update_fragment_mz(assay)
        
        if self.unknown_fragment_type == 'ignore':
            pass
        elif self.unknown_fragment_type == 'keep':
            j = 0
            for i, _ in enumerate(fragments['fragmentMZ']):
                if i in fragment_index:
                    fragments['fragmentMZ'][i] = \
                        assay['fragments']['fragmentMZ'][j]
                    j += 1
            
            assay['fragments'] = fragments
                
        assay = self.assay_builder.filter_fragments_by_mz(
            assay
        )
            
        return assay
        
    
    def decoy(self, assay):
        new_assay = copy.deepcopy(assay)        
        new_assay.update(self.decoy_sequence(new_assay))
        new_assay = self.update_fragments(new_assay)
        
        metadata = new_assay.get('metadata', None)
        if metadata is None:
            new_assay['metadata'] = {'decoy': True}
        else:
            metadata['decoy'] = True        
            if metadata.get('protein', None) is not None:
                metadata['protein'] = '/'.join((
                    'DECOY_' + x
                    for x in str(metadata['protein']).split('/')
                ))
        return new_assay
        
    
        
if __name__ == '__main__':
    decoy = DecoyAssayGenerator(fragment_types=[])
    
    assay0 = {
        'peptideSequence': 'TNVSKFDLPNR', 
        'modification': None, 
        'precursorCharge': 3, 
        'precursorMZ': 998.42809, 
        'glycanStruct': '(N(N(H(H(H(H(H(H(H(H))))))))))', 
        'fragments': {
            'fragmentMZ': [216.097882035, 315.166296035, 402.19832403500004, 530.293287035, 677.361701035, 792.388644035, 905.472708035, 175.11895173499994, 289.1618787349999, 386.2146427349999, 499.29870673499994, 614.325649735, 761.3940637349999, 889.4890267349999, 976.5210547349999, 1075.5894687349999, 299.13499584000004, 398.20340984, 760.39881484, 875.42575784, 988.50982184, 1085.5625858399999, 1272.6695095399998, 445.24815138499997, 488.76416538499996, 538.298372385, 595.319835885, 543.2849309375, 636.8383927875, 1493.759447235, 1696.8388197349998, 1858.8916431349999, 1373.7171885399998, 747.383361635, 848.923047885, 929.949459585, 1010.9758712849999, 1092.002282985, 1173.028694685, 1254.0551063849998, 1335.081518085, 1416.1079297849997, 687.3622322875, 498.59133310166663, 566.2844572683333, 620.3020650683333, 674.3196728683332, 728.3372806683334, 782.3548884683333, 836.3724962683333], 
            'fragmentNumber': [2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 3, 6, 7, 8, 9, 10, 7, 8, 9, 10, 9, 10, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None], 
            'fragmentGlycan': [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, 'Y-N(1)', 'Y-N(2)', 'Y-N(2)H(1)', 'Y$', 'Y-N(1)', 'Y-N(2)', 'Y-N(2)H(1)', 'Y-N(2)H(2)', 'Y-N(2)H(3)', 'Y-N(2)H(4)', 'Y-N(2)H(5)', 'Y-N(2)H(6)', 'Y-N(2)H(7)', 'Y$', 'Y-N(1)', 'Y-N(2)', 'Y-N(2)H(1)', 'Y-N(2)H(2)', 'Y-N(2)H(3)', 'Y-N(2)H(4)', 'Y-N(2)H(5)'], 
            'fragmentLossType': ['noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss', 'noloss'], 
            'fragmentType': ['b', 'b', 'b', 'b', 'b', 'b', 'b', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'y', 'b$', 'b$', 'b$', 'b$', 'b$', 'b$', 'y$', 'y', 'y', 'y', 'y', 'b$', 'y$', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y'], 
            'fragmentAnnotation': ['b2^+1', 'b3^+1', 'b4^+1', 'b5^+1', 'b6^+1', 'b7^+1', 'b8^+1', 'y1^+1', 'y2^+1', 'y3^+1', 'y4^+1', 'y5^+1', 'y6^+1', 'y7^+1', 'y8^+1', 'y9^+1', 'b$2^+1', 'b$3^+1', 'b$6^+1', 'b$7^+1', 'b$8^+1', 'b$9^+1', 'y$10^+1', 'y7^+2', 'y8^+2', 'y9^+2', 'y10^+2', 'b$9^+2', 'y$10^+2', 'Y-N(1)^+1', 'Y-N(2)^+1', 'Y-N(2)H(1)^+1', 'Y$^+1', 'Y-N(1)^+2', 'Y-N(2)^+2', 'Y-N(2)H(1)^+2', 'Y-N(2)H(2)^+2', 'Y-N(2)H(3)^+2', 'Y-N(2)H(4)^+2', 'Y-N(2)H(5)^+2', 'Y-N(2)H(6)^+2', 'Y-N(2)H(7)^+2', 'Y$^+2', 'Y-N(1)^+3', 'Y-N(2)^+3', 'Y-N(2)H(1)^+3', 'Y-N(2)H(2)^+3', 'Y-N(2)H(3)^+3', 'Y-N(2)H(4)^+3', 'Y-N(2)H(5)^+3'], 
            'fragmentIntensity': [1009609.5, 562854.7, 139883.1, 89470.6, 103570.8, 283033.8, 205681.9, 393638.4, 129843.6, 2256995, 5616303, 1028598, 1680987.1, 514180.9, 1157310.5, 577459.3, 146154.4, 102992.3, 27188.2, 83698.5, 32845.3, 42232.4, 70465.8, 24330.3, 170790.8, 153452.1, 73777.4, 21220.6, 101923.3, 8737480, 464471.7, 351479.3, 783311.3, 17168948, 5641570, 5025644, 3937760.5, 3747871.8, 3250004, 3960476.8, 2553953.3, 420593.9, 646249.6, 81078.5, 1249016.5, 1184480.6, 396399.9, 199690.1, 194890.7, 104705.2], 
            'fragmentCharge': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3]
        }
    }

    print(decoy.decoy(assay0))