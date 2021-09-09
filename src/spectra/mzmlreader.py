from pymzml.run import Reader

class MzmlReader():
    def __init__(self, file, parameters=None):
        self.reader = Reader(file)
        
        if parameters is None:
            parameters = mzml_reader_parameters()
        self.parameters = parameters
        
        
    def __enter__(self):
        return self        
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        pass
        
    
    def read_spectrum(self):
        spec = next(self.reader, None)
        if spec is None:
            return None               
        
        def set_value(result, path, value):
            if isinstance(path, str):
                result[path] = value
            elif isinstance(path, list):
                d = result
                for i, field in enumerate(path): 
                    if i == len(path) - 1:
                        d[field] = value
                        break
                    d = d.setdefault(field, {})
        
        def set_data(result, params, spec): 
            mz = list(spec.mz)
            intensity = list(spec.i)
            
            param = params.get('mz', None)
            if param is not None:
                path = param.get('path', None)
                set_value(result, path, mz)
                
            param = params.get('intensity', None)
            if param is not None:
                path = param.get('path', None)
                set_value(result, path, intensity)
                
        
        def set_params(result, params, spec):            
            if params is None:
                return
            
            for parameter in params:
                path = parameter.get('path', None)
                if path is None:
                    continue
                
                value = None
                name = parameter.get('name', None)
                if name is not None:
                     value = spec.get(name, None)
                
                if value is None:
                    accession = parameter.get('accession', None)
                    if accession is not None:
                         cv = spec.ms.get(accession, None)
                         if cv is not None:
                             value = cv.get('value', None)
                 
                if value is not None:                       
                    convert = parameter.get('convert', None)
                    if callable(convert):
                        value = convert(value)
                        
                if value is None:
                    func = parameter.get('function', None)
                    if callable(func):
                        value = func(spec)
                
                if value is not None:
                    set_value(result, path, value)
                    
        if spec.get('total ion current chromatogram', None) is not None or \
            spec.get('MS:1000235', None) is not None:
            return self.read_spectrum()        
        
        result = {}        
        ms1 = spec.get('MS1 spectrum', None) is not None or \
            spec.ms.get('MS:1000579', None) is not None
        
        params = self.parameters.get('params', None)
        if isinstance(params, dict):
            if ms1:
                params = params.get('MS1', None)
            else:
                params = params.get('MSn', None)            
        if isinstance(params, list):
            set_params(result, params, spec)
        
        data_params = self.parameters.get('data', None)
        if isinstance(data_params, dict):
            if ms1:
                data_params = data_params.get('MS1', data_params)
            else:
                data_params = data_params.get('MSn', data_params)            
        if isinstance(data_params, dict):
            set_data(result, data_params, spec)
                
        return result
        
    
def mzml_reader_parameters():
    def spectrum_ref(spec):
        precursor_element = spec.xmlTreeIterFree.find(
            'precursorList/precursor'
        )
        if precursor_element is None:
            precursor_element = spec.xmlTreeIterFree.find(
                '{http://psi.hupo.org/ms/mzml}precursorList' +
                '/{http://psi.hupo.org/ms/mzml}precursor'
            )
        
        if precursor_element is not None:
            spectrum_ref = precursor_element.get('spectrumRef')
        else:
            spectrum_ref = None
        return spectrum_ref
    
    def activation(spec):
        if spec.ms.get('MS:1000133', None) is not None:
            return 'CID'
        elif spec.ms.get('MS:1000422', None) is not None:
            return 'HCD'
        elif spec.ms.get('MS:1000250', None) is not None:
            return 'ECD'
        elif spec.ms.get('MS:1000598', None) is not None:
            return 'ETD'
        elif spec.ms.get('MS:1002631', None) is not None:
            return 'EThcD'
        else:
            return None
    
    parameters = {
        'data': {
            'MS1': {
                'mz': {
                    'path': ['peaks', 'mz']
                },
                'intensity': {
                    'path': ['peaks', 'intensity']
                }
            },
            'MSn': {
                'mz': {
                    'path': ['fragments', 'fragmentMZ']
                },
                'intensity': {
                    'path': ['fragments', 'fragmentIntensity']
                }
            },
        },
        'params': [
            {
                'path': 'msLevel',
                'name': 'ms level',
                'accession': 'MS:1000511',
                'convert': int
            },
            {
                'path': 'rt',
                'name': 'scan start time',
                'accession': 'MS:1000016',
                'convert': float
            },
            {
                'path': ['metadata', 'scanWindowLowerLimit'],
                'name': 'scan window lower limit',
                'accession': 'MS:1000501',
                'convert': float
            },
            {
                'path': ['metadata', 'scanWindowUpperLimit'],
                'name': 'scan window upper limit',
                'accession': 'MS:1000500',
                'convert': float
            },
            {
                'path': ['metadata', 'title'],
                'name': 'spectrum title',
                'accession': 'MS:1000796'
            },
            {
                'path': ['metadata', 'centroid'],
                'name': 'centroid spectrum',
                'accession': 'MS:1000127',
                'convert': lambda x: x is not None
            },
            {
                'path': ['metadata', 'profile'],
                'name': 'profile spectrum',
                'accession': 'MS:1000128',
                'convert': lambda x: x is not None
            },
            {
                'path': ['metadata', 'id'],
                'function': lambda spec: spec.xmlTreeIterFree.get('id', None)
            },
            {
                'path': ['metadata', 'isolationWindowLowerOffset'],
                'name': 'isolation window lower offset',
                'accession': 'MS:1000828',
                'convert': float
            },
            {
                'path': ['metadata', 'isolationWindowTargetMZ'],
                'name': 'isolation window target m/z',
                'accession': 'MS:1000827',
                'convert': float
            },
            {
                'path': ['metadata', 'isolationWindowUpperOffset'],
                'name': 'isolation window upper offset',
                'accession': 'MS:1000829',
                'convert': float
            },
            {
                'path': ['metadata', 'activation'],
                'function': activation
            },
            {
                'path': ['metadata', 'collisionEnergy'],
                'name': 'collision energy',
                'accession': 'MS:1000045',
                'convert': float
            },
            {
                'path': 'precursorMZ',
                'name': 'selected ion m/z',
                'accession': 'MS:1000744',
                'convert': float
            },
            {
                'path': 'precursorIntensity',
                'name': 'peak intensity',
                'accession': 'MS:1000042',
                'convert': float
            },
            {
                'path': 'precursorCharge',
                'name': 'charge state',
                'accession': 'MS:1000041',
                'convert': float
            },        
            {
                'path': ['metadata', 'spectrumRef'],
                'function': spectrum_ref
            },
            {
                'path': ['metadata', 'compensationVoltage'],
                'name': 'FAIMS compensation voltage',
                'accession': 'MS:1001581',
                'convert': float
            }
        ]
    }
    
    return parameters