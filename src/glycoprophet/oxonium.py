import bisect
import os
import click
import itertools
import pandas as pd
from collections import Counter

from pepmass.glycomass import GlycanNode
from assay.oxonium import OxoniumIonExtractor
from spectra.background import BackgroundEstimator
from spectra.mzmlreader import MzmlReader

class OxoniumIonsValidator:
    def __init__(self, data, 
                 background_estimator=True,
                 absolute_intensity=10,
                 relative_intensity=0.05,
                 **kwargs):
        lack_columns = {
            'GlycanStruct', 'filename', 'mz',
            'RT', 'leftWidth', 'rightWidth'
        }.difference(data.columns)
        if len(lack_columns) > 0:
            raise ValueError('missing column(s): ' + str(lack_columns))
        self.data = data
        
        self.extractor = OxoniumIonExtractor(**kwargs)
        
        if background_estimator == True:
            self.background_estimator = BackgroundEstimator()
        elif background_estimator == False:
             self.background_estimator = None
        else:
            self.background_estimator = background_estimator
        self.absolute_intensity = absolute_intensity
        self.relative_intensity = relative_intensity
        
    
    def extract_oxonium_ions(self, spectra):
        oxonium_ions = {}        
        for spec in spectra: 
            if spec.get('msLevel', 2) != 2:
                continue
                 
            oxo = self.extractor.extract_oxonium_ions(spec)
            
            if 'basePeakIntensity' not in oxo['metadata']:
                oxo['metadata']['basePeakIntensity'] = \
                    max(spec['fragments']['fragmentIntensity'])                
            if 'backgroundIntensity' not in oxo['metadata'] and \
                self.background_estimator is not None:
                oxo['metadata']['backgroundIntensity'] = \
                    self.background_estimator.estimate_background(spec)
            
            filename = spec['metadata']['file']
            oxo_dict = oxonium_ions.get(filename, None)
            if oxo_dict is None:
                oxo_dict = oxonium_ions.setdefault(filename, {})
                
            isolation_window = (
                spec['metadata']['isolationWindowTargetMZ'] - \
                    spec['metadata']['isolationWindowLowerOffset'],
                spec['metadata']['isolationWindowTargetMZ'] + \
                    spec['metadata']['isolationWindowUpperOffset']
            )
            oxo_list = oxo_dict.get(isolation_window, None)
            if oxo_list is None:
                oxo_list = oxo_dict.setdefault(isolation_window, [])
            oxo_list.append(oxo)
        
        for oxo_dict in oxonium_ions.values():
            for oxo_list in oxo_dict.values():
                oxo_list.sort(key=lambda x: x['rt'])
        
        self.oxonium_ions = oxonium_ions
    
    
    def validate_oxonium_ions(self):
        rt_indexes = {
            filename: {
                isowin: [x['rt'] for x in oxo_list]
                for isowin, oxo_list in oxo_dict.items()
            }
            for filename, oxo_dict in self.oxonium_ions.items()
        }
        
        def find_nearest_spectrum(filename, precursor_mz, 
                                  rt, rt_window=10, left=None, right=None):
            rt_dict = rt_indexes.get(filename, None)
            if rt_dict is None:
                return None, None
            
            window, rt_list = next(iter(sorted(
                filter(
                    lambda x: x[0][0] <= precursor_mz and x[0][1] >= precursor_mz,
                    rt_dict.items() 
                ), 
                key=lambda x: abs((x[0][0] + x[0][1]) / 2 - precursor_mz)
            )), (None, None))
                        
            if rt_list is None:
                return None, None
            
            index = bisect.bisect_left(rt_list, rt)
            if index == len(rt_list):
                index = index - 1
            elif index == 0:
                pass
            else:
                if abs(rt_list[index] - rt) > abs(rt_list[index - 1] - rt):
                    index = index - 1
            
            if rt_window is not None and abs(rt_list[index] - rt) > rt_window:
                return window, None     
#            if left is not None and rt_list[index] < left:
#                return window, None
#            if right is not None and rt_list[index] > right:
#                return window, None
            return window, index
        
        
        def count_oxonium_ions(oxonium_ions):
            intensity_threshold = max(
                self.absolute_intensity or 0.0,
                oxonium_ions['metadata'].get('basePeakIntensity', 0.0) * \
                    self.relative_intensity \
                    if self.relative_intensity is not None else 0.0,
                oxonium_ions['metadata'].get('backgroundIntensity', 0.0)
            )
            
            count = dict(Counter((
                a.split(':')[0]
                for a, i in zip(
                    oxonium_ions['fragments']['fragmentAnnotation'],
                    oxonium_ions['fragments']['fragmentIntensity']
                )
                if i >= intensity_threshold
            )))            
            count.update({
                k: 0
                for k in set((
                    a.split(':')[0] 
                    for a in self.extractor.oxonium_ions \
                        ['fragments']['fragmentAnnotation']
                )).difference(count.keys())
            })
            
            return count
        
        
        def match_glycan_oxonium_ions(count, glycan_struct):
            for m in GlycanNode.from_str(glycan_struct).composition().keys():
                m_count = [
                    v for k, v in count.items()
                    if m in k.split('+')
                ]
                if len(m_count) > 0 and sum(m_count) == 0:
                    return False
            
            return True
            
        
        def validate_row(row):
            window, index = find_nearest_spectrum(
                filename=row['filename'], precursor_mz=row['mz'],
                rt=row['RT'], left=row['leftWidth'], right=row['rightWidth']
            )
            if window is None or index is None: 
                import warnings
                warnings.warn(
                    'cannot find spectrum filename={filename}, RT={rt}, mz={precursor_mz}' \
                    .format(filename=row['filename'], precursor_mz=row['mz'], rt=row['RT']))
                return pd.Series()
            
            count = count_oxonium_ions(
                oxonium_ions=self.oxonium_ions[row['filename']][window][index]
            )
            matched = match_glycan_oxonium_ions(
                count, glycan_struct=row['GlycanStruct']
            )            
            return pd.Series({'oxonium_matched': matched}) \
                .append(pd.Series(count))
                    
        data = self.data \
            .apply(validate_row, axis=1, result_type='expand') \
            .fillna(-1) \
            .astype(int) 
        self.data = pd.concat((self.data, data), axis=1)        
        return self.data
    
    
def validate_oxonium_ions(report_file, mzml_files, out_file=None, **kwargs):
    if out_file is None:
        out_file = os.path.splitext(report_file)[0] + '_oxonium.tsv'
       
    click.echo("Info: Loading {0}".format(report_file))
    data = pd.read_csv(report_file, sep='\t')
    
    def mzml_as_iterable(path):
        filename = \
            next(filter(
                lambda x: path.endswith(x),
                data['filename'].drop_duplicates()
            ), None) or \
            next(filter(
                lambda x: x.endswith(path),
                data['filename'].drop_duplicates()
            ), None) or \
            path
        
        click.echo("Info: Loading {0}".format(filename))
        
        reader = MzmlReader(path)
        i = 0
        click.echo("Info: Scan", nl=False)  
        while True:
            spec = reader.read_spectrum()        
            if spec is None:   
                break
            spec['metadata']['file'] = filename
            spec['rt'] *= 60
            
            i += 1
            if i % 1000 == 0:
                click.echo(" {0}".format(i), nl=False)                
            yield spec
        
        click.echo("")
        click.echo("Info: {0} loaded, {1} scans".format(filename, i))
    
    
    spectra = itertools.chain.from_iterable((
        mzml_as_iterable(mzml_file)
        for mzml_file in mzml_files
    ))
    
    validator = OxoniumIonsValidator(data=data, **kwargs)
    validator.extract_oxonium_ions(spectra)
    data = validator.validate_oxonium_ions()
    
    data.to_csv(out_file, sep='\t', index=False)
    
    click.echo("Info: Results written {0}".format(out_file))
    
    