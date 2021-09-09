from assay import AssayBuilder

class SpectrumAnnotator():
    def __init__(self, assay_builder=None,
                 tolerance=20, tolerance_unit='ppm', 
                 criteria='mostintense'):
        if assay_builder is None:
            assay_builder = AssayBuilder()
        self.assay_builder = assay_builder
        
        self.tolerance = tolerance
        self.tolerance_unit = tolerance_unit
        self.criteria = criteria

    def annotate(self, spectrum, sequence, modification=None, 
                 **kwargs):
        theoretical = self.assay_builder.theoretical_fragments(
            sequence=sequence, 
            modification=modification,
            **kwargs
        )
        index = match_fragments(
            theoretical, spectrum, 
            tolerance=self.tolerance,
            tolerance_unit=self.tolerance_unit,
            criteria=self.criteria,
            all1=False, all2=True
        )
        
        spectrum['fragments'] = spectrum['fragments'].copy()
        for k, v in theoretical['fragments'].items():
            if k != 'fragmentMZ' and k != 'fragmentIntensity':
                values = spectrum['fragments'].setdefault(k, \
                    [None] * len(spectrum['fragments']['fragmentMZ']))
                for i, j in index:
                    if i is not None and j is not None:
                        values[j] = v[i]

        spectrum.update({
            'peptideSequence': sequence,
            'modification': modification,
            
        })
        return spectrum
            
        
        
def match_peaks(spectrum1, spectrum2, 
                tolerance, tolerance_unit='ppm', 
                criteria='mostintense',
                all1=False, all2=False,
                ms_level=1):    
    if ms_level == 2:
        peaks_key = 'fragments'
        mz_key = 'fragmentMZ'
        intensity_key = 'fragmentIntensity'
    else:
        peaks_key = 'peaks'
        mz_key = 'mz'
        intensity_key = 'intensity'
    
    if tolerance_unit == 'ppm':
        def mz_tolerance_range(mz, tolerance):        
            return (
                mz * (1 - tolerance * 1e-6),
                mz * (1 + tolerance * 1e-6)
            )
    elif tolerance_unit == 'Da':
        def mz_tolerance_range(mz, tolerance):        
            return (
                mz - tolerance,
                mz + tolerance
            )
    else:
        raise ValueError('invalid tolerance unit: ' + str(tolerance_unit))
        
    if criteria == 'mostintense':
        def choose_first(target_index, target_mz, 
                         index1, mz1, index2, mz2):
            intensity = spectrum2[peaks_key].get(intensity_key, None)
            if intensity is not None:
                return intensity[index1] >= intensity[index2]
            else:
                False
    elif criteria == 'nearst':
        def choose_first(target_index, target_mz, 
                         index1, mz1, index2, mz2):
            return abs(mz1 - target_mz) <= abs(mz2 - target_mz)
    else:
        raise ValueError('invalid criteria: ' + str(criteria))
    
    index = []
    used = set()
    for i, x in enumerate(spectrum1[peaks_key][mz_key]):
        l, r = mz_tolerance_range(x, tolerance)        
        j0 = None
        y0 = None
        for j, y in enumerate(spectrum2[peaks_key][mz_key]): 
            if y <= r and y >= l and j not in used:
                if j0 is None or not choose_first(
                    target_index=i, target_mz=x,
                    index1=j0, mz1=y0,
                    index2=j, mz2=y
                ):
                    j0 = j
                    y0 = y
        if j0 is not None:
            index.append((i, j0))
            used.add(j0)
        elif all1:
            index.append((i, None))
        
    if all2:        
        for j, y in enumerate(spectrum2[peaks_key][mz_key]): 
            if j not in used:                
                index.append((None, j))
                
    return index


def match_fragments(spectrum1, spectrum2, 
                  tolerance, tolerance_unit='ppm', 
                  criteria='mostintense',
                  all1=False, all2=False):                    
    return match_peaks(
        spectrum1, spectrum2, 
        tolerance, tolerance_unit=tolerance_unit, 
        criteria=criteria,
        all1=all1, all2=all2,
        ms_level=2
    )
        
            
                
if __name__ == '__main__':
    spec = {'precursorCharge': 3, 'metadata': {'title': 'MouseBrain-Z-T-2.4.4.3.0.dta'}, 'fragments': {'fragmentMZ': [126.05489, 137.33168, 138.05475, 139.05836, 144.06555, 144.81012, 152.22743, 157.64909, 168.06528, 169.0688, 173.38754, 173.42172, 186.07552, 191.17856, 204.08626, 205.0892, 259.88007, 280.05307, 291.02487, 350.31296, 366.13882, 367.14014, 368.14255, 467.69983, 483.60669, 512.19519, 513.20148, 658.25427, 713.27118, 730.30865, 747.22308, 804.31744, 811.82611, 813.34149, 814.33521, 842.35803, 877.88702, 895.88837, 896.90173, 897.88617, 913.89124, 914.89026, 929.84875, 930.84344, 933.37537, 934.38623, 935.36609, 946.41125, 986.37787, 986.8678, 987.38153, 1059.39001, 1067.40955, 1067.90662, 1068.40906, 1068.89209, 1079.43689, 1080.43604, 1140.44226, 1140.93713, 1141.41943, 1141.92432, 1142.4364, 1213.96008, 1214.94397, 1282.52209, 1460.5542, 1461.54456, 1527.2655, 1606.61279, 1607.61975, 1608.62402, 1622.59558, 1623.62061, 1768.66858, 1769.65552, 1770.64307, 1771.67249, 1884.13794, 1914.72058, 1915.73926], 'fragmentIntensity': [5943.0, 388.3, 46186.1, 2890.5, 3310.2, 433.2, 479.3, 498.8, 20506.8, 1073.6, 432.8, 502.8, 8978.8, 570.4, 38373.6, 3577.9, 547.6, 844.4, 570.3, 615.6, 25948.8, 4401.1, 669.0, 764.3, 663.6, 6067.4, 2127.2, 710.8, 918.3, 685.2, 663.9, 738.0, 615.0, 3463.9, 657.4, 719.0, 971.5, 6191.9, 3651.4, 3052.8, 3104.3, 792.8, 823.1, 1872.2, 4462.7, 1171.0, 947.3, 602.0, 2659.5, 1036.1, 1111.8, 2391.2, 4722.1, 5777.4, 3077.7, 994.6, 6050.7, 2506.4, 7508.1, 6945.5, 4154.3, 3459.4, 732.3, 1811.3, 728.9, 679.1, 1205.9, 737.3, 680.4, 4488.7, 3630.2, 921.6, 2279.2, 848.3, 5818.1, 4793.5, 2123.9, 1796.8, 618.9, 751.2, 1187.4]}, 'rt': 1.567459, 'precursorMZ': 931.025978}

    from glycoassay import GlycoAssayBuilder
    annotator = SpectrumAnnotator(
        assay_builder=GlycoAssayBuilder(), tolerance=50
    )
    
    print(annotator.assay_builder.filter_fragments_by_type(annotator.annotate(
        spectrum=spec, 
        sequence='SJMTMWSSQP',
        glycan_struct='(N(N(H(H(N(F)))(H(N(F))))))',
        modification=[
            {'name': 'Oxidation', 'position': 3, 'site': 'M'},
            {'name': 'Oxidation', 'position': 5, 'site': 'M'}
        ]
    )))
    