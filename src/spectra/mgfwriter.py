class MgfWriter():
    def __init__(self, file, parameters=None):
        if isinstance(file, str):
            file = open(file, 'w')
        self.file = file
        
        if parameters is None:
            parameters = mgfwriter_parameters()
        self.parameters = parameters
        
    
    def __enter__(self):
        return self        
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.file.close()
        
    def close(self):
        self.file.close()
        
    
    def write_spectrum(self, spectrum):
        def convert_local_params(spectrum, param):
            value = None            
            path = param.get('path', None)
            if isinstance(path, str):
                if path in spectrum:
                    value = spectrum[path]
            elif isinstance(path, list):
                result = spectrum
                for i, field in enumerate(path): 
                    if i == len(path) - 1 and field in result:
                        value = result[field]
                        break
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
                return func(spectrum)
            
            default = param.get('default', None)            
            return default
        
        lines = ['BEGIN IONS']
        
        local_params = self.parameters.get('localParams', None) 
        if local_params is not None:
            for key, param in local_params.items():
                value = convert_local_params(spectrum, param)
                if value is not None:
                    lines.append(key + '=' + str(value))
        
        ions_params = self.parameters.get('ions', None)
        if ions_params is not None:
            mz_param = ions_params.get('mz')
            mz_array = convert_local_params(spectrum, mz_param)
            intensity_param = ions_params.get('intensity')
            intensity_array = convert_local_params(spectrum, intensity_param)
            
            charge_param = ions_params.get('charge', None)
            if charge_param is not None:
                charge_array = convert_local_params(spectrum, charge_param)
            else:
                charge_array = None
                
            annotation_param = ions_params.get('annotation', None)
            if annotation_param is not None:
                annotation_array = convert_local_params(spectrum, annotation_param)
            else:
                annotation_array = None
                
            for i, mz in enumerate(mz_array):
                line = str(mz) + ' ' + str(intensity_array[i])
                if charge_array is not None:
                    line += ' ' + str(charge_array[i]) + \
                        ('+' if charge_array[i] >= 0 else '')
                if annotation_array is not None:
                    line += ' ' + annotation_array[i]
                lines.append(line)
                
        lines.append('END IONS')
        lines.append('')
        
        self.file.writelines((x + '\n' for x in lines))


def mgfwriter_parameters():
    def pepmass(spectrum):        
        mz = spectrum.get('precursorMZ', None)
        if mz is None:
            return None
        
        line = str(mz)
        
        intensity = spectrum.get('precursorIntensity', None)
        if intensity is not None:
            line += ' ' + str(intensity)
            
        return line
      
    def ion_mobility(spectrum):
        drift_time = spectrum.get('driftTime', None)
        if drift_time is None:
            return None
        
        return pepmass(spectrum) + ' ' + str(drift_time)
    
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
                'convert': lambda x: str(x) + ('+' if x >= 0 else '')
            },
            'RTINSECONDS': {
                'path': 'rt'
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
                'path': ['metadata', 'tolerance']
            },
            'TOLU': {
                'path': ['metadata', 'toleranceUnit']
            }
        }
    }
    
    return parameters



if __name__ == '__main__':
    spec = {
        'fragments': {
            'fragmentMZ': [240.1, 242.1, 245.2, 1623.7, 1624.7], 
            'fragmentIntensity': [3.0, 12.0, 32.0, 55.0, 23.0]
        }, 
        'precursorIntensity': 25674.3, 
        'metadata': {
            'scans': '3', 
            'tolerance': 3.0, 
            'sequence': 'n-AC[DHK]', 
            'title': 'Spectrum 1', 
            'toleranceUnit': 'Da', 
            'composition': '2[H]0[M]3[DE]*[K]'
        }, 
        'precursorMZ': 896.05, 
        'precursorCharge': 3, 
        'rt': 25.0
    }
    
    from io import StringIO
    s = StringIO()
    with MgfWriter(s) as mgf:
        mgf.write_spectrum(spec)
        print(s.getvalue())
        
        