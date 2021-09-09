import numpy as np

class BackgroundEstimator:
    def __init__(self, 
                 noise_peak_ratio=2.0,
                 isotope_number=4, 
                 isotope_tolerance=0.01, 
                 isotope_tolerance_unit='Da'):
        self.noise_peak_ratio = noise_peak_ratio
        self.isotope_number = isotope_number
        self.isotope_tolerance = isotope_tolerance
        self.isotope_tolerance_unit = isotope_tolerance_unit
        
    
    def estimate_background(self, spectrum):
        peak_dict = spectrum.get('fragments', None)
        if peak_dict is not None:
            peaks = np.column_stack((
                peak_dict['fragmentMZ'],
                peak_dict['fragmentIntensity']
            ))
        if peak_dict is None:
            peak_dict = spectrum.get('peaks', None)
            if peak_dict is not None:
                peaks = np.column_stack((
                    peak_dict['mz'],
                    peak_dict['intensity']
                ))     
        if peak_dict is None:
            raise ValueError('missing peaks')
                
        if peaks.shape[0] < 0:
            return 0.0
        
        lower, upper = np.percentile(
            peaks[:, 1], q=(0, 70), 
            interpolation='lower'
        )
        
        if upper <= lower + 0.001:
            return 0.0
                      
        for bk in np.linspace(lower, upper, 21):              
            filtered_peaks = peaks[peaks[:, 1] > bk, :]
            mz_dist = filtered_peaks[1:, 0] - filtered_peaks[:-1, 0]
            int_diff = filtered_peaks[1:, 1] - filtered_peaks[:-1, 1]
            
            if self.isotope_tolerance_unit in {'Da' ,'Th'}:
                isotope_tolerance = self.isotope_tolerance
            elif self.isotope_tolerance_unit == 'ppm':
                isotope_tolerance = filtered_peaks[1:, 0] * \
                    self.isotope_tolerance * 1e-6
            else:
                raise ValueError('invalid isotope_tolerance_unit: ' + \
                                 str(self.isotope_tolerance_unit))
                
            count = [
                sum(
                    (mz_dist > 1.0 / i - isotope_tolerance) &
                    (mz_dist < 1.0 / i + isotope_tolerance) &
                    (int_diff < 0)
                )
                for i in range(1, 1 + self.isotope_number)
            ]
            noise = sum(mz_dist < 1.0 / self.isotope_number - \
                        isotope_tolerance)
            
            if noise < sum(count) * self.noise_peak_ratio:
                return float(bk)
        
        return 0.0
            
            