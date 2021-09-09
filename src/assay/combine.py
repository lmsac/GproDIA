import itertools
import numpy as np

from assay.modseq import stringify_modification
from pepmass.glycomass import GlycanNode

class AssayCombiner():
    def __init__(self, group_key=None):
        if group_key is None:
            group_key = glycopeptide_group_key()
        self.group_key = group_key
        
        
    def combine(self, *assays, return_generator=False):
        return self.remove_redundant(
            itertools.chain.from_iterable(assays),
            return_generator=return_generator
        )
        
    
    def remove_redundant(self, assays, return_generator=False):
        result = (
            self.combine_replicates(spectra)
            for spectra in self.group_replicates(assays)
        )
        result = (x for x in result if x is not None)
        
        if not return_generator:
            result = list(result)
            
        return result
    
        
    def group_replicates(self, assays):
        def get_key(assay):
            return tuple(
                str(k(assay) if callable(k) else assay.get(k, None))
                for k in self.group_key
            )
        
        return (
            list(v)
            for k, v in itertools.groupby(
                sorted(assays, key=get_key), 
                key=get_key
            )
        )
        
    
    def combine_replicates(self, spectra):
        if len(spectra) == 0:
            return None 
        
        return spectra[0]
    
    
class BestReplicateAssayCombiner(AssayCombiner):    
    def __init__(self, 
                 group_key=None,
                 score='score', higher_score_better=True):
        super(BestReplicateAssayCombiner, self) \
            .__init__(group_key=group_key)
        self.score = score
        self.higher_score_better = higher_score_better
        
    
    def combine_replicates(self, spectra):
        if len(spectra) == 0:
            return None 
        
        score = [
            spec['metadata'][self.score]
            for spec in spectra
        ]
        if self.higher_score_better:
            index = np.argmax(score)
        else:
            index = np.argmin(score)
        
        return spectra[index]
    

def glycopeptide_group_key(use_glycan_struct=True, use_glycan_site=True, 
                           within_run=False):
    if use_glycan_struct:
        glycan_key = 'glycanStruct'
    else:
        def glycan_key(x): 
            x = x.get('glycanStruct', None)
            return x and GlycanNode \
                .from_str(x) \
                .composition_str()
    
    group_key = [
        'peptideSequence',
        lambda x: stringify_modification(x.get('modification', None)),
        glycan_key,
        'precursorCharge'
    ]
    
    if use_glycan_site:
        group_key.append('glycanSite')
    
    if within_run:
        def filename(x):
            metadata = x.get('metadata', None)
            if metadata is not None:
                return metadata.get('file', None)
            return None                
        group_key.insert(0, filename)
    
    return group_key
    