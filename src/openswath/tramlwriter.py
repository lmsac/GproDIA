from lxml import etree
from assay.modseq import ModifiedSequenceConverter


class TramlWriter():
    def __init__(self, parameters=None, cv_map=None):
        if parameters is None:
            parameters = traml_writer_parameters()
        self.parameters = parameters
        
        if cv_map is None:
            cv_map = cv_name_accession_map()
        self.cv_map = cv_map
        

    def initialize_document(self):
        self.xmlns_uris = {
            None: 'http://psi.hupo.org/ms/traml',
            'xsi': 'http://www.w3.org/2001/XMLSchema-instance' 
        }
        
        self.root = etree.Element(
            'TraML',
            attrib={ 
                etree.QName(self.xmlns_uris['xsi'], 'schemaLocation'): \
                    'http://psi.hupo.org/ms/traml TraML1.0.0.xsd'
            },
            version='1.0.0',
            nsmap=self.xmlns_uris
        )
        
        self.cv_list = etree.SubElement(self.root, 'cvList')
        etree.SubElement(
            self.cv_list, 'cv', 
            attrib={ 
                'id': 'MS',
                'fullName': 'Proteomics Standards Initiative Mass Spectrometry Ontology',
                'version': 'unknown',
                'URI': 'http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo'
            }
        )
        etree.SubElement(
            self.cv_list, 'cv', 
            attrib={ 
                'id': 'UO',
                'fullName': 'Unit Ontology',
                'version': 'unknown',
                'URI': 'http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo'
            }
        )
        
        self.protein_list = etree.SubElement(self.root, 'ProteinList')
        self.compound_list = etree.SubElement(self.root, 'CompoundList')
        self.transition_list = etree.SubElement(self.root, 'TransitionList')
        
        self.protein_dict = {}
        self.peptide_dict = {}
        self.transition_dict = {}
        return self.root
        
    
    def create_param(self, name, value=None, unit_name=None, **attribs):
        name = str(name)
        accession = self.cv_map.get(name, None)
        if unit_name is not None:
            unit_name = str(unit_name)
            unit_accession = self.cv_map.get(unit_name, None)
            if unit_accession is None:
                raise ValueError('unknown unit: ' + str(unit_name))
            attribs.update({
                'unitCvRef': unit_accession.split(':', 1)[0],
                'unitAccession': unit_accession,
                'unitName': unit_name
            })
           
        if value is not None:
            attribs.update({
                'value': str(value)
            })
            
        if accession is not None:
            attribs.update({
                'cvRef': accession.split(':', 1)[0],
                'accession': accession,
                'name': name
            })            
            return etree.Element('cvParam', attrib=attribs)
        else:
            attribs.update({                
                'name': name
            })
            if 'type' not in attribs:
                attribs.update({
                    'type': 'xsd:string'
                })            
            return etree.Element('userParam', attrib=attribs)
        
    
    def add_protein(self, id, protein_accession, sequence, params=None):
        id = str(id)
        if id in self.protein_dict:
            raise ValueError('protein already exists: ' + id)
        
        protein = etree.Element('Protein', attrib={'id': id})
        
        if protein_accession is not None:
            param = {
                'name': 'protein accession',
                'value': protein_accession
            }
            if params is None:
                params = [param]
            else:
                params = [param] + [
                    x for x in params 
                    if x.get('name', None) != 'protein accession'
                ]
            
        
        if params is not None:
            for x in params:
                protein.append(
                    self.create_param(**x)
                )
        
        seq = etree.SubElement(protein, 'Sequence')
        if sequence is not None:
            seq.text = str(sequence)
        
        self.protein_list.append(protein)
        self.protein_dict[id] = protein
        return protein
    
    
    def add_peptide(self, id, sequence, protein_ref, 
                    charge_state=None,
                    params=None, 
                    rt=None, rt_unit=None, rt_list=None):        
        id = str(id)
        if id in self.peptide_dict:
            raise ValueError('peptide already exists: ' + id)
            
        peptide = etree.Element('Peptide', attrib={
            'id': id, 
            'sequence': str(sequence)
        })
        
        if charge_state is not None:            
            param = {
                'name': 'charge state',
                'value': charge_state
            }
            if params is None:
                params = [param]
            else:
                params = [param] + [
                    x for x in params 
                    if x.get('name', None) != 'charge state'
                ]
            
        if params is not None:
            for x in params:
                peptide.append(
                    self.create_param(**x)
                )
        
        etree.SubElement(peptide, 'ProteinRef', attrib={'ref': str(protein_ref)})
        
        if rt is not None:
            rt_list = [[{
                'name': 'normalized retention time', 
                'value': rt,
                'unit_name': rt_unit
            }]]
                    
        if rt_list is not None:
            retention_time_list = etree.SubElement(peptide, 'RetentionTimeList')
            for rt_params in rt_list:
                retention_time = etree.SubElement(retention_time_list, 
                                                  'RetentionTime')
                for x in rt_params:
                    retention_time.append(
                        self.create_param(**x)
                    )
           
        self.compound_list.append(peptide)
        self.peptide_dict[id] = peptide
        return peptide
    
    
    def add_transition(self, id, peptide_ref, 
                       precursor_mz=None,
                       product_mz=None, product_intensity=None,
                       product_charge=None,
                       precursor_params=None, product_params=None, 
                       interpretation_list=None,
                       params=None):
        id = str(id)
        if id in self.transition_dict:
            raise ValueError('transition already exists: ' + id)
            
        transition = etree.Element('Transition', attrib={
            'id': id, 
            'peptideRef': str(peptide_ref)
        })
                
        if precursor_mz is not None:
            param = {
                'name': 'isolation window target m/z',
                'value': precursor_mz, 
                'unit_name': 'm/z'
            }
            if precursor_params is None:
                precursor_params = [param]
            else:
                precursor_params = [param] + [
                    x for x in precursor_params 
                    if x.get('name', None) != 'isolation window target m/z'
                ]
                
        if product_intensity is not None:
            param = {
                'name': 'product ion intensity',
                'value': product_intensity
            }
            if params is None:
                params = [param]
            else:
                params = [param] + [
                    x for x in params 
                    if x.get('name', None) != 'product ion intensity'
                ]
                
        if product_mz is not None:
            param = {
                'name': 'isolation window target m/z',
                'value': product_mz
            }
            
            if product_params is None:
                product_params = [param]
            else:
                product_params = [param] + [
                    x for x in product_params 
                    if x.get('name', None) != 'isolation window target m/z'
                ]   
            
        if product_charge is not None:
            param = {
                'name': 'charge state',
                'value': product_charge
            }
            
            if product_params is None:
                product_params = [param]
            else:
                product_params = [param] + [
                    x for x in product_params 
                    if x.get('name', None) != 'charge state'
                ]
    
        if precursor_params is not None:
            precursor = etree.SubElement(transition, 'Precursor')
            for x in precursor_params:
                precursor.append(
                    self.create_param(**x)
                )
            
        product = etree.SubElement(transition, 'Product')
        if product_params is not None:            
            for x in product_params:
                product.append(
                    self.create_param(**x)
                )
        
        if interpretation_list is not None:
            interp_list = etree.SubElement(product, 'InterpretationList')
            for interp_params in interpretation_list:
                interp = etree.SubElement(interp_list, 
                                          'Interpretation')
                for x in interp_params:
                    interp.append(
                        self.create_param(**x)
                    )
                    
        if params is not None:
            for x in params:
                transition.append(
                    self.create_param(**x)
                )
        
        self.transition_list.append(transition)
        self.transition_dict[id] = transition
        return transition
    

    def add_assay(self, assay, **kwargs):
        def get_value(param, assay, **kwargs):
            if param is None:
                return None
            value = None
            path = param.get('path', None)
            if isinstance(path, str):
                if path in assay:
                    return assay[path]
            elif isinstance(path, list):
                result = assay
                for i, field in enumerate(path): 
                    if i == len(path) - 1 and field in result:
                        return result[field]
                    result = result.get(field, None)
                    if result is None:
                        break
            
            if value is not None:            
                convert = param.get('convert', None)
                if callable(convert):
                    value = convert(value)
                return value
            
            func = param.get('function', None)
            if callable(func):
                return func(assay, **kwargs)
            
            default = param.get('default', None)            
            return default
        
        def get_params(params, assay, **kwargs):
            if params is None:
                return None
            
            def get_param(param, assay, **kwargs):
                return {
                    k: (get_value(v, assay, **kwargs) \
                        if isinstance(v, dict) \
                        else v)
                    for k, v in param.items()
                }
               
            return list(filter(
                lambda x: None not in x.values(), 
                (
                    get_param(param, assay, **kwargs)    
                    for param in params
                )
            ))
                
        def get_transition_params(transition_params, i):
            return list(filter(
                lambda x: None not in x.values(), 
                ({
                    k: (v[i] if isinstance(v, list) else v)
                    for k, v in p.items()
                }
                for p in transition_params)
            ))
                             
            
        protein = self.parameters.get('protein', None)
        if protein is not None:
            protein_id = get_value(protein.get('id'), assay, **kwargs)
            
            if protein_id not in self.protein_dict:
                protein_args = {
                    k: get_value(v, assay, **kwargs)
                    for k, v in protein.items()
                    if k != 'params' and k != 'id'
                }
                protein_params = get_params(
                    protein.get('params', None), assay, 
                    **kwargs
                )
                
                self.add_protein(
                    id=protein_id, 
                    params=protein_params, 
                    **protein_args
                )        
        else:
            protein_id = None            
           
        peptide = self.parameters.get('peptide', None)
        if peptide is not None:
            peptide_id = get_value(peptide.get('id'), assay, **kwargs)
            
            if peptide_id not in self.peptide_dict:
                peptide_args = {
                    k: get_value(v, assay, **kwargs)
                    for k, v in peptide.items()
                    if k != 'params' and k != 'id'
                }
                peptide_params = get_params(
                    peptide.get('params', None), assay, 
                    **kwargs
                )
                
                self.add_peptide(
                    id=peptide_id,
                    protein_ref=protein_id,
                    params=peptide_params,
                    **peptide_args
                )                
        else:
            peptide_id = None            
            
        transition = self.parameters.get('transition', None)
        if transition is not None:
            transition_id = get_value(transition.get('id'), assay, **kwargs)
            
            transition_args = {
                k: get_value(v, assay, **kwargs)
                for k, v in transition.items()
                if k != 'params' and k != 'id' and k != 'interpretation'
            }
            transition_params = get_params(
                transition.get('params', None), assay, 
                **kwargs
            )
            interpretation = get_params(
                transition.get('interpretation', None), assay, 
                **kwargs
            )
            
            for i, x in enumerate(transition_id):
                self.add_transition(
                    id=x,
                    peptide_ref=peptide_id,
                    params=get_transition_params(transition_params, i),
                    interpretation_list=[
                        get_transition_params(interpretation, i)
                    ],
                    **get_transition_params([transition_args], i)[0]
                )
       
        
    def assays_to_traml(self, assays, **kwargs):
        self.initialize_document()
        
        for i, x in enumerate(assays):
            self.add_assay(x, index=i, **kwargs)
        
        return self.root
    
    
    def save_traml(self, file, pretty_print=True, **kwargs):
        s = etree.tostring(self.root, pretty_print=pretty_print, **kwargs)
        
        if isinstance(file, str):
            file = open(file, 'wb')
        
        file.write(s)
        
        file.close()
        
        
        
def cv_name_accession_map():
    return {
        'protein accession': 'MS:1000885',
        'm/z': 'MS:1000040',
        'charge state': 'MS:1000041',
        'isolation window target m/z': 'MS:1000827',
        'peptide group label': 'MS:1000893',
        
        'normalized retention time': 'MS:1000896',
        'iRT retention time normalization standard': 'MS:1002005',         
        
        'product ion series ordinal': 'MS:1000903',        
        'product interpretation rank': 'MS:1000926',
        'product ion intensity': 'MS:1001226',
        'product ion m/z delta': 'MS:1000904',        
            
        'frag: y ion'           : 'MS:1001220',
        'frag: b ion - H2O'     : 'MS:1001222',
        'frag: y ion - H2O'     : 'MS:1001223',
        'frag: b ion'           : 'MS:1001224',
        'frag: x ion'           : 'MS:1001228',
        'frag: a ion'           : 'MS:1001229',
        'frag: z ion'           : 'MS:1001230',
        'frag: c ion'           : 'MS:1001231',
        'frag: b ion - NH3'     : 'MS:1001232',
        'frag: y ion - NH3'     : 'MS:1001233',
        'frag: a ion - H2O'     : 'MS:1001234',
        'frag: a ion - NH3'     : 'MS:1001235',
        'frag: d ion'           : 'MS:1001236',
        'frag: v ion'           : 'MS:1001237',
        'frag: w ion'           : 'MS:1001238',
        'frag: immonium ion'    : 'MS:1001239',
        'frag: internal yb ion' : 'MS:1001365',
        'frag: internal ya ion' : 'MS:1001366',
        'frag: z+1 ion'         : 'MS:1001367',
        'frag: z+2 ion'         : 'MS:1001368',
        'frag: c ion - H2O'     : 'MS:1001515',
        'frag: c ion - NH3'     : 'MS:1001516',
        'frag: z ion - H2O'     : 'MS:1001517',
        'frag: z ion - NH3'     : 'MS:1001518',
        'frag: x ion - H2O'     : 'MS:1001519',
        'frag: x ion - NH3'     : 'MS:1001520',
        'non-identified ion'    : 'MS:1001240',
                
        'target SRM transition': 'MS:1002007',
        'decoy SRM transition': 'MS:1002008',
        
        'second': 'UO:0000010',
        'minute': 'UO:0000031',        
    }


def traml_writer_parameters(): 
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
                }
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
                    'name': 'annotation', 
                    'value': {
                        'path': ['fragments', 'fragmentAnnotation']
                    }
                }
            ],
        }
    }
            