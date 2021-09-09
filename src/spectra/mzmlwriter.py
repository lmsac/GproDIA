import psims.mzml
import numpy as np

class MzmlWriter:
    def __init__(self, file, parameters=None):
        self.writer = psims.mzml.MzMLWriter(file)
        
        if parameters is None:
            parameters = mzml_writer_parameters()
        self.parameters = parameters
        
        
    def __enter__(self):
        self.writer.__enter__()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.spectrum_list.__exit__(exc_type, exc_val, exc_tb)
        self.run.__exit__(exc_type, exc_val, exc_tb)
        self.writer.__exit__(exc_type, exc_val, exc_tb)        
    
    
    def begin(self):        
        self.writer.begin()
        self.writer.controlled_vocabularies()
        
        self.writer.file_description([
            'MS1 spectrum',
            'MSn spectrum'
        ])
        self.writer.software_list([{
            'id': 'psims-writer',
            'version': psims.__version__,
        }])
    
        self.writer.instrument_configuration_list([
            self.writer.InstrumentConfiguration(
                id="IC1", 
                component_list=[
                    self.writer.Source(1, [
                        'electrospray ionization', 
                        'electrospray inlet'
                    ]), 
                    self.writer.Analyzer(2, [
                        'orbitrap'
                    ]),
                    self.writer.Detector(3, [
                        'inductive detector'
                    ])
                ]
            )
        ])
                    
        self.writer.data_processing_list([
            self.writer.DataProcessing([
                self.writer.ProcessingMethod(
                    order=1, software_reference='psims-writer', 
                    params=['']
                )
            ], id='DP1')
        ])
                
        self.run = self.writer.run(id=1, instrument_configuration='IC1')
        self.run.__enter__()
        
        self.spectrum_list = self.writer.spectrum_list(count=65535)
        self.spectrum_list.__enter__()
        self.spectrum_count = 0
        
    
    def close(self):
        self.spectrum_list.__exit__(None, None, None)
        self.run.__exit__(None, None, None)
        self.writer.end()
        self.writer.close()
        
        
    def write_spectrum(self, spectrum):
        def get_value(param, spectrum):
            if param is None:
                return None
            value = None
            path = param.get('path', None)
            if isinstance(path, str):
                if path in spectrum:
                    return spectrum[path]
            elif isinstance(path, list):
                result = spectrum
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
                return func(spectrum)
            
            default = param.get('default', None)            
            return default
        
        params = self.parameters.get('params', None)
        data_params = self.parameters.get('data', None)
        
        ms_level = spectrum.get('msLevel', None)
        if ms_level == 1:
            if isinstance(params, dict):
                params = params.get('MS1', params)
            if isinstance(data_params, dict):
                data_params = data_params.get('MS1', data_params) 
        else:
            if isinstance(params, dict):
                params = params.get('MSn', params)
            if isinstance(data_params, dict):
                data_params = data_params.get('MSn', data_params)
                        
        spec_args = {}
        if ms_level is not None:
            spec_args.setdefault('params', [{'ms level': ms_level}])        
        
        if data_params is not None:
            for key, param in data_params.items():
                value = get_value(param, spectrum)
                if value is None:
                    continue
                
                if key in {'mz', 'intensity', 'charge'}:
                    spec_args[key + '_array'] = \
                        np.array(value)
                
        if params is not None:
            for key, param in params.items():
                value = get_value(param, spectrum)
                if value is None:
                    continue
                                
                if key in {'id', 'polarity', 'centroided', 'scan_start_time'}: 
                    spec_args[key] = value
                
                elif key == 'scan_window_lower_limit': 
                    spec_args.setdefault('scan_window_list', [[None, None]]) \
                        [0][0] = value
                elif key == 'scan_window_upper_limit':
                    spec_args.setdefault('scan_window_list', [[None, None]]) \
                        [0][1] = value
                
                elif key == 'spectrum_ref':
                    spec_args.setdefault('precursor_information', {}) \
                        ['spectrum_reference'] = value
                elif key == 'selected_ion_mz':
                    spec_args.setdefault('precursor_information', {}) \
                        ['mz'] = value
                elif key == 'peak_intensity':
                    spec_args.setdefault('precursor_information', {}) \
                        ['intensity'] = value
                elif key == 'charge_state':
                    spec_args.setdefault('precursor_information', {}) \
                        ['charge'] = value
                elif key == 'isolation_window_lower_offset':
                    spec_args.setdefault('precursor_information', {}) \
                        .setdefault('isolation_window_args', {}) \
                        ['lower'] = value
                elif key == 'isolation_window_target_mz':
                    spec_args.setdefault('precursor_information', {}) \
                        .setdefault('isolation_window_args', {}) \
                        ['target'] = value
                elif key == 'isolation_window_upper_offset':
                    spec_args.setdefault('precursor_information', {}) \
                        .setdefault('isolation_window_args', {}) \
                        ['upper'] = value
                elif key == 'activation':
                    spec_args.setdefault('precursor_information', {}) \
                        .setdefault('activation', []) \
                        .append(value)     
                elif key == 'collision_energy':
                    spec_args.setdefault('precursor_information', {}) \
                        .setdefault('activation', []) \
                        .append({'collision energy': value}) 
                else:
                    spec_args.setdefault('params', []) \
                        .append({key: value})
        
        intensity_array = spec_args.get('intensity_array', None)
        if intensity_array is not None:
            spec_args.setdefault('params', [])  \
                .append({'total ion current': intensity_array.sum()})
                
        spec_args.setdefault('id', self.spectrum_count)
                
        self.writer.write_spectrum(**spec_args)
        self.spectrum_count += 1
        
        
def mzml_writer_parameters(): 
    parameters = {
        'data': {
            'MS1': {
                'mz': {
                    'path': ['peaks', 'mz']
                },
                'intensity': {
                    'path': ['peaks', 'intensity']
                },
                'charge': {
                    'path': ['peaks', 'charge']
                }
            },
            'MSn': {
                'mz': {
                    'path': ['fragments', 'fragmentMZ']
                },
                'intensity': {
                    'path': ['fragments', 'fragmentIntensity']
                },
                'charge': {
                    'path': ['fragments', 'fragmentCharge']
                }
            },
        },
        'params': {
            'id': {
                'path': ['metadata', 'id']
            },
            'polarity': {
                'path': ['metadata', 'polarity'],                
            },
            'centroided': {
                'path': ['metadata', 'centroid']
            },
            'scan_start_time': {
                'path': 'rt'
            },
            'scan_window_lower_limit': {
                'path': ['metadata', 'scanWindowLowerLimit']
            },
            'scan_window_upper_limit': {
                'path': ['metadata', 'scanWindowUpperLimit']
            },
            'spectrum_ref': {
                'path': ['metadata', 'spectrumRef'],
            },
            'selected_ion_mz': {
                'path': 'precursorMZ',
            },
            'peak_intensity': {
                'path': 'precursorIntensity'
            },
            'charge_state': {
                'path': 'precursorCharge'
            },
            'activation': {
                'path': ['metadata', 'activation']
            },
            'collision_energy': {
                'path': ['metadata', 'collisionEnergy']
            },
            'isolation_window_lower_offset': {
                'path': ['metadata', 'isolationWindowLowerOffset']
            },
            'isolation_window_target_mz': {
                'path': ['metadata', 'isolationWindowTargetMZ']
            },
            'isolation_window_upper_offset': {
                'path': ['metadata', 'isolationWindowUpperOffset']
            },
            'FAIMS compensation voltage': {
                'path': ['metadata', 'compensationVoltage']
            }
        }
    }
    
    return parameters
   
    
    