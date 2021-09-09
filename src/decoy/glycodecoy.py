import copy
import random

from assay import GlycoAssayBuilder
from decoy import DecoyAssayGenerator

class GlycoDecoyAssayGenerator(DecoyAssayGenerator):
    def __init__(self, assay_builder=None, 
                 **kwargs):        
        if assay_builder is None:
            assay_builder = GlycoAssayBuilder()
        self.assay_builder = assay_builder
        
        super(GlycoDecoyAssayGenerator, self) \
            .__init__(assay_builder=assay_builder, **kwargs)
         
    def decoy_sequence(self, assay):
        glycan_site = assay.get('glycanSite', None)
        if glycan_site is not None:
            assay = assay.copy()
            modification = assay.get('modification', None)
            if modification is None:
                modification = []
            assay['modification'] = modification + [{
                'name': 'GlycoMod', 'position': glycan_site
            }]
    
        result = super(GlycoDecoyAssayGenerator, self) \
            .decoy_sequence(assay)
        
        if glycan_site is not None:
            t = next((
                (i, x['position'])
                for i, x in enumerate(result['modification'])
                if x['name'] == 'GlycoMod'
            ))
            result['glycanSite'] = t[1]
            result['modification'].pop(t[0])
        
        return result
            
        
    def peptide_decoy(self, assay):        
        assay = super(GlycoDecoyAssayGenerator, self) \
            .decoy(assay=assay)
        
        metadata = assay.get('metadata', None)
        if metadata is None:
            assay['metadata'] = {'decoy': True, 'peptideDecoy': True}
        else:
            metadata['peptideDecoy'] = True 
        return assay
        
        
    def glycan_decoy(self, assay):
        new_assay = copy.deepcopy(assay)
        
        if self.unknown_fragment_type == 'ignore':
            new_assay = self.assay_builder.filter_fragments_by_type(new_assay)
            
        fragments = new_assay['fragments']
                
        for i, x in enumerate(fragments['fragmentType']):
            if x in self.assay_builder.fragment_types:
                continue
            elif x in self.assay_builder.glycan_fragment_types:   
                if fragments['fragmentGlycan'][i] == x + '0' or \
                    fragments['fragmentGlycan'][i] == x + '$':
                    continue
                    
                mz_shift = round(random.uniform(1, 30), 2) / \
                    (fragments['fragmentCharge'][i] or 1)                
                fragments['fragmentMZ'][i] += mz_shift                 
                fragments['fragmentType'][i] = 'DECOY_' + x                
                fragments['fragmentAnnotation'][i] = 'DECOY_' + \
                    fragments['fragmentAnnotation'][i] + '[+' + \
                    '{:.2f}'.format(mz_shift) + ']'
            elif self.unknown_fragment_type == 'error':
                raise ValueError('fragment not found: ' + x)
                
#        precursor_mz = new_assay.get('precursorMZ', None)
#        if precursor_mz is not None:
#            new_assay['precursorMZ'] = precursor_mz + \
#                round(random.uniform(1, 30), 2) / \
#                new_assay.get('precursorCharge', 1)
        
        metadata = new_assay.get('metadata', None)
        if metadata is None:
            new_assay['metadata'] = {'decoy': True, 'glycanDecoy': True}
        else:
            metadata['decoy'] = True  
            metadata['glycanDecoy'] = True  
        return new_assay
    
    
    def decoy(self, assay, decoy_type='peptide'):
        if decoy_type == 'peptide':
            return self.peptide_decoy(assay)
        elif decoy_type == 'glycan':
            return self.glycan_decoy(assay)
        
    
        
if __name__ == '__main__':
    decoy = GlycoDecoyAssayGenerator()
    
    assay0 = {
        'peptideSequence': 'TJVSKFDLPNR', 
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
    
    print(decoy.decoy(assay0, decoy_type='glycan'))
    
    