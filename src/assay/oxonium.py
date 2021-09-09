import copy
from pepmass.glycomass import GlycoPeptideMassCalculator
from assay.annotation import match_fragments


class OxoniumIonExtractor:
    def __init__(self, 
                 mass_calculator=None,
                 mz_tolerance=20, mz_tolerance_unit='ppm',
                 **kwargs):
        if mass_calculator is None:
            mass_calculator = GlycoPeptideMassCalculator()
        self.mass_calculator = mass_calculator

        self.mz_tolerance = mz_tolerance
        self.mz_tolerance_unit = mz_tolerance_unit
        
        self.calculate_oxonium_ion_mz()
        
    
    def calculate_oxonium_ion_mz(self):
        oxonium_ion_mz = self.mass_calculator.oxonium_ion_mz()
        self.oxonium_ions = {
            'fragments': {
                'fragmentMZ': oxonium_ion_mz['fragment_mz'],
                'fragmentAnnotation': oxonium_ion_mz['fragment_name']
            }
        }


    def extract_oxonium_ions(self, spectrum):
        index = match_fragments(
            spectrum, self.oxonium_ions, 
            tolerance=self.mz_tolerance, tolerance_unit=self.mz_tolerance_unit
        )
        spectrum = copy.deepcopy(spectrum)   
        spectrum.update({    
            'fragments': {
                'fragmentMZ': \
                    [spectrum['fragments']['fragmentMZ'][i] \
                     for i, j in index],
                'fragmentIntensity': \
                    [spectrum['fragments']['fragmentIntensity'][i] \
                     for i, j in index],
                'fragmentAnnotation': \
                    [self.oxonium_ions['fragments']['fragmentAnnotation'][j] \
                     for i, j in index]
                     
            }
        })
        return spectrum
    
        
        