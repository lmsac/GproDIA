class MgfReader():
    def __init__(self, file, parameters=None):
        if isinstance(file, str):
            file = open(file, 'r')
        self.file = file
        
        if parameters is None:
            parameters = mgfreader_parameters()
        self.parameters = parameters
        
    
    def __enter__(self):
        return self        
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()
        
    def close(self):
        self.file.close()
        
    
    def read_spectrum(self):
        local_params = self.parameters.get('localParams', None)        
        ions_params = self.parameters.get('ions', None)
        
        def set_local_params(spectrum, key, value):
            parameter = None
            if local_params is not None:
                parameter = local_params.get(key, None)

            if parameter is None:
                from warnings import warn
                warn('unknown parmeter: ' + key)
                return
                
            path = parameter.get('path', None)
            if path is not None:
                convert = parameter.get('convert', None)
                if callable(convert):
                    value = convert(value)
                
                if isinstance(path, str):
                    spectrum[path] = value
                elif isinstance(path, list):
                    d = spectrum
                    for i, field in enumerate(path): 
                        if i == len(path) - 1:
                            d[field] = value
                            break
                        d = d.setdefault(field, {})
                return
            
            func = parameter.get('function', None)
            if callable(func):
                func(spectrum, value)
            
        def set_ions(spectrum, mz, intensity, charge=None, annotation=None):
            if ions_params is None:
                return
            
            def get_array_or_create(key):
                parameter = ions_params.get(key, None)
                if parameter is not None:
                    path = parameter.get('path', None)
                    if isinstance(path, str):
                        l = spectrum.get(path, None)
                        if l is None:
                            l = spectrum.setdefault(path, [])
                        return l
                    elif isinstance(path, list):
                        d = spectrum
                        for i, field in enumerate(path): 
                            if i == len(path) - 1:
                                l = d.get(field, None)
                                if l is None:
                                    l = d.setdefault(field, [])
                                break
                            d = d.setdefault(field, {})
                        return l
            
            mz_array = get_array_or_create('mz')
            if mz_array is not None:
                mz_array.append(mz)
            intensity_array = get_array_or_create('intensity')
            if intensity_array is not None:
                intensity_array.append(intensity)
            if charge is not None:
                charge_array = get_array_or_create('charge')
                if charge_array is not None:
                    charge_array.extend([None] * \
                        (len(mz_array) - len(charge_array)))
                    charge_array[-1] = charge
            if annotation is not None:
                annotation_array = get_array_or_create('annotation')
                if annotation_array is not None:
                    annotation_array.extend([None] * \
                        (len(mz_array) - len(annotation_array)))
                    annotation_array[-1] = annotation
        
        spectrum = None
        while True:
            line = self.file.readline()
            if not line:
                if spectrum is not None:
                    raise ValueError('unexpected EOF')
                return None
            
            line = line.strip()
            if len(line) == 0 or \
                line[0] in {'#', ';', '!', '/'}:
                continue
            
            if line == 'BEGIN IONS':
                if spectrum is not None:
                    raise ValueError('[Offset ' + str(self.file.tell()) + \
                        '] invalid format: ' + line)
                spectrum = {}
                continue
            
            if line == 'END IONS':
                if spectrum is None:
                    raise ValueError('[Offset ' + str(self.file.tell()) + \
                        '] invalid format: ' + line)
                return spectrum
            
            if spectrum is None:
                continue
            
            s = line.split('=', 1)
            if len(s) == 2:                                
                set_local_params(spectrum, key=s[0], value=s[1])
                continue
                    
            s = line.split(' ', 4) 
            if len(s) < 2:
                raise ValueError('[Offset ' + str(self.file.tell()) + \
                    '] invalid format: ' + line)
            if len(s) >= 3:
                charge = int(s[2].strip('+'))  
            else:
                charge = None
            if len(s) == 4:
                annotation = s[3]
            else:
                annotation = None
           
            set_ions(
                spectrum, 
                mz=float(s[0]), 
                intensity=float(s[1]), 
                charge=charge, 
                annotation=annotation
            )
            
                         
def mgfreader_parameters():
    def pepmass(spectrum, value):
        s = value.split(' ')
        if len(s) == 0 or len(s) > 3:
            raise ValueError('invalid format: ' + 'PEPMASS=' + value)

        spectrum['precursorMZ'] = float(s[0])
        
        intensity = None
        charge = None
        if len(s) == 2:
            if s[1].endswith('+'):
                charge = s[1]
            else:
                intensity = s[1]
        if len(s) == 3:
            intensity = s[1]
            charge = s[2]

        if intensity is not None:
            spectrum['precursorIntensity'] = float(intensity)
        if charge is not None:
            spectrum.setdefault('precursorCharge', int(charge.strip('+')))
      
    def ion_mobility(spectrum, value):
        s = value.split(' ')
        if len(s) <= 1:
            raise ValueError('invalid format: ' + 'ION_MOBILITY=' + value)

        spectrum['driftTime'] = float(s[-1])
    
    parameters = {
        'ions': {
            'mz': {
                'path': ['fragments', 'fragmentMZ']
            },
            'intensity': {
                'path': ['fragments', 'fragmentIntensity']
            },
            'charge': {
                'path': ['fragments', 'fragmentCharge']
            },
            'annotation': {
                'path': ['fragments', 'fragmentAnnotation']
            },
        },
        'localParams': {
            'TITLE': {
                'path': ['metadata', 'title']
            },
            'PEPMASS': {
                'function': pepmass
            },
            'CHARGE': {
                'path': 'precursorCharge',
                'convert': lambda x: int(x.strip('+'))
            },
            'RTINSECONDS': {
                'path': 'rt',
                'convert': float
            },
            'ION_MOBILITY': {
                'function': ion_mobility
            },
            
            'COMP': {
                'path': ['metadata', 'composition' ]
            },
            'ETAG': {
                'path': ['metadata', 'errorTolerantSequenceTag' ]
            },
            'INSTRUMENT': {
                'path': ['metadata', 'instrument']
            },
            'IT_MODS': {
                'path': ['metadata', 'variableModifications']
            },
            'LOCUS': {
                'path': ['metadata', 'locus']
            },
            'RAWFILE': {
                'path': ['metadata', 'rawFile']
            },
            'RAWSCANS': {
                'path': ['metadata', 'rawScans']
            },
            'SCANS': {
                'path': ['metadata', 'scans']
            },
            'SEQ': {
                'path': ['metadata', 'sequence']
            },
            'TAG': {
                'path': ['metadata', 'sequenceTag']
            },
            'TOL': {
                'path': ['metadata', 'tolerance'],
                'convert': float
            },
            'TOLU': {
                'path': ['metadata', 'toleranceUnit']
            }
        }
    }
    
    return parameters
    

    
if __name__ == '__main__':
    s = '''BEGIN IONS
TITLE=Spectrum 1
PEPMASS=896.05 25674.3
SCANS=3
RTINSECONDS=25
CHARGE=3+
TOL=3
TOLU=Da
SEQ=n-AC[DHK]
COMP=2[H]0[M]3[DE]*[K]
240.1 3
242.1 12
245.2 32
1623.7 55
1624.7 23
END IONS
'''
    
    from io import StringIO
    with MgfReader(StringIO(s)) as mgf:
        while True:
            spec = mgf.read_spectrum()
            if spec is not None:
                print(spec)
            else:
                break
    
    