import copy
import re

from .assay import AssayBuilder
from pepmass import GlycoPeptideMassCalculator

class GlycoAssayBuilder(AssayBuilder):
    def __init__(self, 
                 mass_calculator=None, 
                 glycan_fragment_types=None,
                 glycan_fragment_charges=None,
                 **kwargs):
        if mass_calculator is None:
            mass_calculator = GlycoPeptideMassCalculator()
        
        fragment_types = kwargs.pop('fragment_types', None)
        if fragment_types is None:
            fragment_types = ['b', 'y', 'b-N(1)', 'y-N(1)', 'b$', 'y$']
        fragment_loss_types = kwargs.pop('fragment_loss_types', None)
        if fragment_loss_types is None:
            fragment_loss_types = ['noloss']

        super(GlycoAssayBuilder, self).__init__(
            mass_calculator=mass_calculator,
            fragment_types=fragment_types,
            fragment_loss_types=fragment_loss_types,            
            **kwargs
        )
        
        if glycan_fragment_types is None:
            glycan_fragment_types = ['Y']
        self.glycan_fragment_types = glycan_fragment_types
        
        if glycan_fragment_charges is None:
            glycan_fragment_charges = [1, 2, 3]
        self.glycan_fragment_charges = glycan_fragment_charges          
        
        
    def assay(self, sequence, charge, glycan_struct=None, glycan_site=None,
              modification=None, fragments=None, 
              **kwargs):
        result = super(GlycoAssayBuilder, self).assay(
            sequence=sequence,
            charge=charge,
            modification=modification,
            fragments=fragments,
            **kwargs
        )
        result.update({
            'glycanStruct': glycan_struct,
            'glycanSite': glycan_site
        })
        return result
        
    
    def theoretical_fragments(self, sequence, modification=None, 
                             glycan_struct=None, glycan_site=None,
                             **kwargs):
        def get_fragment_annotation(fragment_type, charge, 
                                fragment_glycan=None, **kwargs):
            loss = kwargs.get('loss', None)        
            return fragment_glycan + \
                ('-' + loss if loss is not None and loss != 'noloss' else '') + \
                '^+' + str(charge)
        
        frag_type = kwargs.pop('fragment_type', None)
        if frag_type is None:
            pep_frag_type = list(self.fragment_types)
            gly_frag_type = list(self.glycan_fragment_types)
        else:
            if isinstance(frag_type, str):
                frag_type = [frag_type]
            pep_frag_type = [
                x for x in frag_type 
                if x in self.fragment_types
            ]
            gly_frag_type = [
                x for x in frag_type 
                if x in self.glycan_fragment_types
            ]        
        
        if glycan_struct is None or \
            isinstance(glycan_struct, str) and \
            not glycan_struct.startswith('(N'):            
            pep_frag_type = [
                x for x in pep_frag_type
                if not (x.endswith('-N(1)') or x.endswith('$'))
            ]                   
        
        if len(pep_frag_type) == 0:
            result = super(GlycoAssayBuilder, self) \
                .assay(
                     sequence=sequence,
                     modification=modification                     
                )        
        else:
            result = super(GlycoAssayBuilder, self) \
                .theoretical_fragments(
                     sequence=sequence,
                     modification=modification,
                     glycan=glycan_struct,
                     glycan_site=glycan_site,
                     fragment_type=pep_frag_type,
                     **kwargs
                )        
           
        if len(gly_frag_type) == 0:
            return result
            
        if glycan_struct is None:
            return result
            
        fragment_mz = []
        fragment_type = []
        fragment_number = []
        fragment_charge = []
        fragment_loss_type = []
        fragment_annotation = []
        fragment_glycan = []

        glycan_fragments = self.mass_calculator.fragment_mz(
            sequence=sequence,
            modification=modification,
            glycan=glycan_struct,
            glycan_site=glycan_site,
            fragment_type=gly_frag_type,
            charge=self.glycan_fragment_charges
        )
                
        for x in glycan_fragments:
            frgtype = x['fragment_type']            
            frgcharge = x['charge']
                        
            frgmz = x['fragment_mz']
            frgglycan = x['fragment_name']
            
            frgannot = [
                get_fragment_annotation(
                    fragment_type=frgtype,
                    fragment_glycan=x,
                    charge=frgcharge
                )
                for x in frgglycan
            ]

            fragment_mz.extend(frgmz)
            fragment_number.extend([None] * len(frgglycan))
            fragment_type.extend([frgtype] * len(frgglycan))
            fragment_charge.extend([frgcharge] * len(frgglycan))
            fragment_loss_type.extend(['noloss'] * len(frgglycan))
            fragment_glycan.extend(frgglycan)
            fragment_annotation.extend(frgannot)
            
        result.update({            
            'glycanStruct': glycan_struct,
            'glycanSite': glycan_site
        })
        result['fragments']['fragmentMZ'].extend(fragment_mz)
        result['fragments']['fragmentType'].extend(fragment_type)
        result['fragments']['fragmentNumber'].extend(fragment_number)
        result['fragments']['fragmentCharge'].extend(fragment_charge)
        result['fragments']['fragmentLossType'].extend(fragment_loss_type)
        result['fragments']['fragmentAnnotation'].extend(fragment_annotation)
        result['fragments']['fragmentGlycan'] = \
            [None] * (len(result['fragments']['fragmentNumber']) - \
            len(fragment_glycan)) + \
            fragment_glycan
            
        return result
        
     
    def update_precursor_mz(self, assay):
        sequence = assay['peptideSequence']
        modification = assay.get('modification', None)
        glycan_struct = assay.get('glycanStruct', None)
        precursor_charge = int(assay['precursorCharge'])
        
        precursor_mz = self.mass_calculator.precursor_mz(
            sequence=sequence,
            modification=modification,
            glycan=glycan_struct,
            charge=precursor_charge
        )
        
        assay.update({
            'precursorMZ': precursor_mz
        })
        return assay
    

    def parse_glycan_fragment_name(self, fragment_name):        
        for gfrg in self.glycan_fragment_types:
            mat = re.match('^[' + gfrg + '][0$]$', fragment_name)
            if mat is not None:
                return fragment_name, None
            
            mat = re.match('^([' + gfrg + '])-(([A-Za-z]+\\([0-9]+\\))+)$', 
                           fragment_name)
            if mat is not None:
                return mat.group(1), { 
                    t[0]:int(t[1]) 
                    for t in re.findall('([A-Za-z]+)\\(([0-9]+)\\)', mat.group(2)) 
                }
            
            mat = re.match('^([' + gfrg + '])-(([A-Za-z]+[0-9]+)+)$', 
                           fragment_name)
            if mat is not None:
                return mat.group(1), { 
                    t[0]:int(t[1]) 
                    for t in re.findall('([A-Za-z]+)([0-9]+)', mat.group(2)) 
                }
                
        raise ValueError('invalid fragment name:' + str(fragment_name))
        
        
    def update_fragment_mz(self, assay, **kwargs):
        fragments = copy.deepcopy(assay['fragments'])
        
        fragment_index = self.filter_fragments_by_type(
            assay, self.glycan_fragment_types, return_index=True
        )
        
        assay = self.filter_fragments_by_index(
            assay, fragment_index=fragment_index, invert=True
        )
        
        glycan_struct = assay.get('glycanStruct', None)              
        glycan_site = assay.get('glycanSite', None)
           
        assay = super(GlycoAssayBuilder, self).update_fragment_mz(
            assay=assay,
            glycan=glycan_struct,
            glycan_site=glycan_site,
            **kwargs
        )
        
        if len(fragment_index) > 0:
            sequence = assay['peptideSequence']
            modification = assay.get('modification', None)            
            
            fragment_type = fragments['fragmentType']
            fragment_glycan = fragments['fragmentGlycan']
            fragment_charge = fragments['fragmentCharge']
            fragment_loss_type = fragments['fragmentLossType']
            
            glycan_fragments = self.mass_calculator.fragment_mz(
                sequence=sequence,
                glycan=glycan_struct,
                glycan_site=glycan_site,
                modification=modification,
                fragment_type=list(set(fragment_type[i] \
                    for i in fragment_index)),
                loss=list(set(fragment_loss_type[i] for i in fragment_index)),
                charge=list(set(fragment_charge[i] for i in fragment_index)),
                **kwargs
            )

            fragment_mz = []
            j = 0                
            for i, _ in enumerate(fragment_type):
                if i not in fragment_index:
                    fragment_mz.append(assay['fragments']['fragmentMZ'][j])
                    j += 1
                    continue
                      
                mz = None
                for x in glycan_fragments:
                    if x['fragment_type'] == fragment_type[i] and \
                        x['charge'] == fragment_charge[i] and \
                        (x.get('loss', 'noloss') == fragment_loss_type[i] or \
                         x.get('loss', 'noloss') == 'noloss' and \
                         fragment_loss_type[i] is None):
                        glycan_composition = self \
                            .parse_glycan_fragment_name(fragment_glycan[i])
                        index = next((
                            i for i, y in enumerate(x['fragment_name'])
                            if self.parse_glycan_fragment_name(y) == \
                                glycan_composition
                        ), None)  
                        if index is None:
                            mz = None
                        else:
                            mz = x['fragment_mz'][index]
                        fragment_mz.append(mz)
                        if mz is None:
                            mz = 0
                        break
                if mz is None:
                    raise ValueError('fragment not found: ', fragment_type[i])
         
            fragments.update({
                'fragmentMZ': fragment_mz
            })
            assay['fragments'] = fragments

        return assay
        

    def filter_fragments_by_type(self, assay, fragment_type=None, 
                                 return_index=False):        
        if fragment_type == None:
            fragment_type = self.fragment_types + self.glycan_fragment_types
        
        return super(GlycoAssayBuilder, self).filter_fragments_by_type(
            assay, 
            fragment_type=fragment_type,
            return_index=return_index
        )
        
        
    def filter_fragments_by_amino_acid_number(self, assay, 
                                              min_amino_acid_number, 
                                              return_index=False):
        fragment_index = super(GlycoAssayBuilder, self) \
            .filter_fragments_by_amino_acid_number(
                assay,
                min_amino_acid_number=min_amino_acid_number,
                return_index=True
            )
            
        fragment_index_1 = self \
            .filter_fragments_by_type(
                assay,
                fragment_type=self.glycan_fragment_types,
                return_index=True
            )
        
        fragment_index = list(sorted(set(fragment_index) \
                                     .union(fragment_index_1)))
        
        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay, 
                fragment_index=fragment_index
            )
      
        
    def filter_glycan_fragments_by_monosaccharide_number(
        self, assay, 
        min_monosaccharide_number, 
        return_index=False
    ):
        def enough_monosaccharide_number(fragment_glycan):
            if fragment_glycan is None:
                monosaccharide_number = 0
            else:
                monosaccharide_number = \
                    self.parse_glycan_fragment_name(fragment_glycan)[1]
                if monosaccharide_number is None:
                    monosaccharide_number = 0
            
            if isinstance(min_monosaccharide_number, dict):
                if isinstance(monosaccharide_number, dict):
                    for k, v in min_monosaccharide_number.items():
                        if k == 'all':
                            if sum(monosaccharide_number.values()) < v:
                                return False
                        elif monosaccharide_number.get(k, 0) < v:
                            return False
                    return True
                else:
                    for k, v in min_monosaccharide_number.items():
                        if monosaccharide_number < v:
                            return False
                    return True
            else:
                if isinstance(monosaccharide_number, dict):
                    return sum(monosaccharide_number.values()) >= \
                        min_monosaccharide_number
                else:
                    return monosaccharide_number >= min_monosaccharide_number
        
        fragment_index = [
            i 
            for i, x in enumerate(assay['fragments']['fragmentGlycan'])
            if enough_monosaccharide_number(x)
        ]
        
        fragment_index_1 = self \
            .filter_fragments_by_type(
                assay,
                fragment_type=self.glycan_fragment_types,
                return_index=True
            )
        
        fragment_index = list(sorted( \
            set(range(0, len(assay['fragments']['fragmentType']))) \
                .difference(fragment_index_1) \
                .union(fragment_index)
        ))
        
        if return_index:
            return fragment_index
        else:
            return self.filter_fragments_by_index(
                assay, 
                fragment_index=fragment_index
            )
        
            
    def filter_fragments(self, assay, 
                         prior_peptide_fragment_number=None,
                         prior_peptide_fragment_criteria=None,
                         prior_glycan_fragment_number=None,
                         prior_glycan_fragment_criteria=None,
                         min_fragment_monosaccharide_number=None,
                         max_fragment_number=None,
                         min_relative_fragment_intensity=None,
                         return_index=False,
                         **kwargs):        
        fragment_index_0 = super(GlycoAssayBuilder, self) \
            .filter_fragments(assay, return_index=True, **kwargs) 
        
        if min_fragment_monosaccharide_number is not None:
            fragment_index_1 = self \
                .filter_glycan_fragments_by_monosaccharide_number(
                    assay, return_index=True,
                    min_monosaccharide_number=min_fragment_monosaccharide_number
                ) 
            fragment_index_0 = list(sorted(set(fragment_index_0) \
                                     .intersection(fragment_index_1)))
        
        assay = super(GlycoAssayBuilder, self) \
            .filter_fragments_by_index(assay, fragment_index_0) 
        
        if prior_peptide_fragment_number is None and \
            prior_glycan_fragment_number is None:
            fragment_index = super(GlycoAssayBuilder, self) \
                .filter_fragments(
                    assay,                     
                    max_fragment_number=max_fragment_number,
                    min_relative_fragment_intensity= \
                        min_relative_fragment_intensity,
                    return_index=True
                )   
                
            if return_index:
                return [fragment_index_0[i] for i in fragment_index]
            else:
                assay = super(GlycoAssayBuilder, self) \
                    .filter_fragments_by_index(
                        assay, 
                        fragment_index
                    )
                return assay
        
        if prior_peptide_fragment_number is not None:
            pep_frag_index = set(self.filter_fragments_by_type(
                assay, fragment_type=self.fragment_types,
                return_index=True
            ))  
            if prior_peptide_fragment_criteria is not None and \
                len(prior_peptide_fragment_criteria) > 0:
                pep_frag_index = pep_frag_index.intersection(
                    self.filter_fragments(
                        assay, **prior_peptide_fragment_criteria,
                        return_index=True
                    )
                )
            
        if prior_glycan_fragment_number is not None:
            gly_frag_index = set(self.filter_fragments_by_type(
                assay, fragment_type=self.glycan_fragment_types,
                return_index=True
            ))
            if prior_glycan_fragment_criteria is not None and \
                len(prior_glycan_fragment_criteria) > 0:
                gly_frag_index = gly_frag_index.intersection(
                    self.filter_fragments(
                        assay, **prior_glycan_fragment_criteria,
                        return_index=True
                    )
                )                
           
        if max_fragment_number is not None:
            fragment_intensity = assay['fragments']['fragmentIntensity']
            fragment_index = sorted(
                range(len(fragment_intensity)), 
                key=lambda k: fragment_intensity[k],
                reverse=True
            )            
            
            pep_frag_index_1 = []
            gly_frag_index_1 = []
            other_frag_index_1 = []            
            for i, x in enumerate(fragment_index):
                if prior_peptide_fragment_number is not None and \
                    x in pep_frag_index and \
                    len(pep_frag_index_1) < prior_peptide_fragment_number:
                    pep_frag_index_1.append(x)
                    
                elif prior_glycan_fragment_number is not None and \
                    x in gly_frag_index and \
                    len(gly_frag_index_1) < prior_glycan_fragment_number:
                    gly_frag_index_1.append(x)
                    
                else:
                    other_frag_index_1.append(x)
            
            other_frag_num = max_fragment_number - len(pep_frag_index_1) - \
                len(gly_frag_index_1)
            if other_frag_num >= 0:
                fragment_index = pep_frag_index_1 + gly_frag_index_1 + \
                    other_frag_index_1[:other_frag_num]
            else:
                pep_frag_num = int(round(max_fragment_number * \
                    prior_peptide_fragment_number / \
                    (prior_peptide_fragment_number + \
                     prior_glycan_fragment_number)))
                gly_frag_num = max_fragment_number - pep_frag_num            
                fragment_index = pep_frag_index_1[:pep_frag_num] + \
                    gly_frag_index_1[:gly_frag_num]
            
        if fragment_index is not None:
            assay_1 = self.filter_fragments_by_index(
                assay, 
                fragment_index=fragment_index
            )            
        else:
            assay_1 = assay
            
        if min_relative_fragment_intensity is not None:
            fragment_index_1 = super(GlycoAssayBuilder, self) \
                .filter_fragments(
                    assay_1, 
                    min_relative_fragment_intensity= \
                        min_relative_fragment_intensity,
                    return_index=True
                )
                 
            if fragment_index is None:
                fragment_index = fragment_index_1
            else:
                fragment_index = [
                    fragment_index[i] 
                    for i in fragment_index_1
                ]
        
        if return_index:
            if fragment_index is None:
                return fragment_index_0
            else:
                return [fragment_index_0[i] for i in fragment_index]
        else:
            if fragment_index is not None:
                assay = super(GlycoAssayBuilder, self) \
                    .filter_fragments_by_index(
                        assay, 
                        fragment_index
                    )
            return assay
        
        
    def filter_assay(self, assay,
                     min_peptide_fragment_number=None,
                     min_peptide_fragment_criteria=None,
                     min_glycan_fragment_number=None,
                     min_glycan_fragment_criteria=None,
                     **kwargs):
        kwargs.setdefault(
            'prior_peptide_fragment_number', 
            min_peptide_fragment_number
        )
        kwargs.setdefault(
            'prior_peptide_fragment_criteria', 
            min_peptide_fragment_criteria
        )
        kwargs.setdefault(
            'prior_glycan_fragment_number',
            min_glycan_fragment_number
        )
        kwargs.setdefault(
            'prior_glycan_fragment_criteria',
            min_glycan_fragment_criteria
        )
                    
        assay = super(GlycoAssayBuilder, self) \
            .filter_assay(assay, **kwargs) 
        
        if assay is None:
            return None
        
        if min_peptide_fragment_number is not None:
            pep_frag_index = set(self.filter_fragments_by_type(
                assay, fragment_type=self.fragment_types,
                return_index=True
            ))
            if min_peptide_fragment_criteria is not None:
                pep_frag_index = pep_frag_index.intersection(
                    self.filter_fragments(
                        assay, 
                        **min_peptide_fragment_criteria,
                        return_index=True
                    )
                )                
            if len(pep_frag_index) < min_peptide_fragment_number:
                return None
        
        if min_glycan_fragment_number is not None:
            gly_frag_index = set(self.filter_fragments_by_type(
                assay, fragment_type=self.glycan_fragment_types,
                return_index=True
            ))
            if min_glycan_fragment_criteria is not None:                
                gly_frag_index = gly_frag_index.intersection(
                    self.filter_fragments(
                        assay, 
                        **min_glycan_fragment_criteria,
                        return_index=True
                    )
                )                
            if len(gly_frag_index) < min_glycan_fragment_number:
                return None         
            
        return assay 
    
        
if __name__ == '__main__':
    assay = GlycoAssayBuilder()
    
    print(assay.theoretical_fragments(
        sequence='TJVSKFDLPNR',
        glycan_struct='(N(N(H(H(H(H(H(H(H(H))))))))))'
    ))