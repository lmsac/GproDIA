from pepmass import ModifiedPeptideMassCalculator
import re
        
class ModifiedSequenceConverter:
    def __init__(self, mass_calculator=None):
        if mass_calculator is not None:
            self.mass_calculator = mass_calculator
        else:
            self.mass_calculator = ModifiedPeptideMassCalculator()
            
    def to_tpp_format(self, sequence, modification=None):
        aa = ['n'] + list(sequence) + ['c']        
        if modification is not None:
            for modsite in modification:
                mod = self.mass_calculator.find_var_mod(sequence, modsite)
                if mod is not None:
                    if modsite.get('site', None) == 'N-term':
                        aa[0] += '['+ str(round(mod.delta_mass + \
                            self.mass_calculator.element_mass('H'))) + ']'
                    elif modsite.get('site', None) == 'C-term':
                        aa[-1] = '['+ str(round(mod.delta_mass + \
                            self.mass_calculator.element_mass('O') + \
                            self.mass_calculator.element_mass('H'))) + ']'
                    else:
                        aa[modsite['position']] += '[' + str(round( \
                            mod.delta_mass + \
                            self.mass_calculator.aa_residue_mass(
                                sequence[modsite['position'] - 1]
                            ))) + ']'
        
        if self.mass_calculator.fixed_modifications is not None:
            for mod in self.mass_calculator.fixed_modifications:
                if mod.site == 'N-term':
                    aa[0] += '['+ str(round(mod.delta_mass + \
                            self.mass_calculator.element_mass('H'))) + ']'
                elif mod.site == 'C-term':
                    aa[-1] = '['+ str(round(mod.delta_mass + \
                        self.mass_calculator.element_mass('O') + \
                        self.mass_calculator.element_mass('H'))) + ']'
                else:
                    for i, x in enumerate(aa):
                        if mod.site == x or \
                            (isinstance(mod.site, list) and x in mod.site):
                            aa[i] += '[' + str(round( \
                                mod.delta_mass + \
                                self.mass_calculator.aa_residue_mass(
                                    sequence[i - 1]
                                ))) + ']'
        
        if aa[0] == 'n':
            aa[0] = ''
        if aa[-1] == 'c':
            aa[-1] = ''
        return ''.join(aa)
    
    
    def parse_tpp_format(self, modified_sequence, include_fixed=True):
        def find_mod(site, mass):
            if mass is None or mass == '':
                return None
            
            if isinstance(mass, str) and \
                (mass.startswith('+') or mass.startswith('-')):
                delta_mass = float(mat[2])
            elif isinstance(mass, float) and mass < 0:
                delta_mass = mass
            else:
                if site == 'N-term':
                    unmod_mass = self.mass_calculator.element_mass('H')
                elif  site == 'C-term':
                    unmod_mass = self.mass_calculator.element_mass('O') + \
                        self.mass_calculator.element_mass('H')
                else:
                    unmod_mass = self.mass_calculator.aa_residue_mass(site)
                
                delta_mass = float(mass) - unmod_mass
                    
            mods = []
            if self.mass_calculator.variable_modifications is not None:
                mods.extend(self.mass_calculator.variable_modifications)
            if include_fixed and \
                self.mass_calculator.fixed_modifications is not None:
                mods.extend(self.mass_calculator.fixed_modifications)
                
            best_mod = None
            best_diff = None
            for mod in mods:
                if mod.site == site or \
                    (isinstance(mod.site, list) and site in mod.site): 
                    diff = mod.delta_mass - delta_mass
                    if abs(diff) <= 0.5 and \
                        (best_diff is None or abs(diff) < abs(best_diff)):
                        best_diff = diff
                        best_mod = mod
                      
            if best_mod is None:
                raise ValueError('modification not found: ' + \
                                 site + '[' + str(mass) + ']')
            return {
                'name': best_mod.name,
                'site': site
            }                
                    
        matches = re.findall(
            '(([A-Za-z])(?:\\[([+\\-]?[0-9]*(?:\\.[0-9]*)?)\\])?)', 
            modified_sequence
        )
        
        aa = []
        n_mods = []
        c_mods = []
        aa_mods = []
        has_n = matches[0][1] == 'n'
        has_c = matches[-1][1] == 'c'
        if has_n:
            if matches[0][2] != '':
                mod = find_mod(site='N-term', mass=matches[0][2])
                mod['position'] = None
                n_mods.append(mod)
            matches.pop(0)
        if has_c:
            if matches[0][2] != '':
                mod = find_mod(site='C-term', mass=matches[-1][2])
                mod['position'] = None
                c_mods.append(mod)
            matches.pop(-1)
        
        for i, mat in enumerate(matches):
            aa.append(mat[1])
            if mat[2] != '':                
                mod = find_mod(site=mat[1], mass=mat[2])
                mod['position'] = i + 1
                aa_mods.append(mod)                
        
        sequence = ''.join(aa)
        modification = n_mods + aa_mods + c_mods
        return sequence, modification
        

def stringify_modification(modification):                
    if modification is None or modification == 'null' or \
        modification == 'None':
        return None
        
    def _to_str(m):
        r = m['name'] + '('
        if m['position'] is not None:
            r += str(m['position'])
        if m['site'] is not None:
            r += str(m['site'])
        r += ')'
        return r
                    
    return ';'.join ((
        _to_str(m)
        for m in modification
    ))
    
    
def parse_modification(s):         
    if s is None or \
        s == 'null' or \
        s == 'None':
        return None
    
    def _from_str(s):                    
        mat = re.match('^(.+)\\(([0-9]+)?([A-Z]|N-term|C-term)?\\)$', s)
        name = mat.group(1)
        pos = mat.group(2)
        if pos is not None:
            pos = int(pos)
        site = mat.group(3)
        return {
            'name': name, 
            'position': pos, 
            'site': site
        }
        
    return [
        _from_str(x)
        for x in s.split(';')
        if x != ''
    ]
    
    
        
if __name__ == '__main__':
    sequence_converter = ModifiedSequenceConverter()
    
    print(sequence_converter.to_tpp_format(
        sequence='PEPCPEPMPEPR',
        modification=[
            {'name': 'Acetyl', 'site': 'N-term'},
            {'name': 'Oxidation', 'position': 8 }
        ]
    ))
    
    print(sequence_converter.parse_tpp_format(
        modified_sequence='n[43]PEPC[160]PEPM[147]PEPR'
    ))