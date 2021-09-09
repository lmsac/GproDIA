import re

class pGlycoToAssayConverter:
    def __init__(self, matched_ion_has_mz=False):
        self.matched_ion_has_mz = matched_ion_has_mz
    
    def parse_modification(self, modification):
        if not isinstance(modification, str) or modification == 'null':
            return None
        
        def _parse_modification(s):
            mat = re.match('^([0-9]+),([^\\[]+)\\[([A-Za-z0-9\\-]+)\\]$', s)
            if mat is None:
                raise ValueError('invalid modification: ' + s)
            pos = int(mat.group(1))
            name = mat.group(2)
            site = mat.group(3)
            if site == 'N-term' or site == 'ProteinN-term' or \
                site == 'PeptideN-term' or site == 'AnyN-term':
                pos = None
                site = 'N-term'
            elif site == 'C-term' or site == 'ProteinC-term' or \
                site == 'PeptideC-term' or site == 'AnyC-term':
                pos = None
                site = 'C-term'
            return {
                'name': name, 
                'position': pos, 
                'site': site
            }
        
        return [
            _parse_modification(s)
            for s in modification.split(';') 
            if s != ''
        ]
    
    
    def parse_ion_annotation(self, label):
        mat = re.match('^([by])([0-9]+)\\+([0-9]+)$', label)
        if mat is not None:
            return {
                'type': 'peptide',
                'fragment_type': mat.group(1), 
                'fragment_number': int(mat.group(2)), 
                'loss': 'noloss',
                'charge': int(mat.group(3))
            }
        
        mat = re.match('^([by])(-N|-N\\(1\\))([0-9]+)\\+([0-9]+)$', label)
        if mat is not None:
            return {
                'group': 'peptide',
                'fragment_type': mat.group(1) + '-N(1)', 
                'fragment_number': int(mat.group(3)), 
                'loss': 'noloss',
                'charge': int(mat.group(4))
            }            
        mat = re.match('^([by])([0-9]+)(-N|-N1|-N\\(1\\))\\+([0-9]+)$', \
                      label)
        if mat is not None:
            return {
                'group': 'peptide',
                'fragment_type': mat.group(1) + '-N(1)', 
                'fragment_number': int(mat.group(2)), 
                'loss': 'noloss',
                'charge': int(mat.group(4))
            }
        
        mat = re.match('^([by])(\\$)([0-9]+)\\+([0-9]+)$', label)
        if mat is not None:
            return {
                'group': 'peptide',
                'fragment_type': mat.group(1) + mat.group(2), 
                'fragment_number': int(mat.group(3)), 
                'loss': 'noloss',
                'charge': int(mat.group(4))
            }            
        mat = re.match('^([by])([0-9]+)(\\$)\\+([0-9]+)$', \
                      label)
        if mat is not None:
            return {
                'group': 'peptide',
                'fragment_type': mat.group(1) + mat.group(3), 
                'fragment_number': int(mat.group(2)), 
                'loss': 'noloss',
                'charge': int(mat.group(4))
            }
        
        mat = re.match('^([Y])\\+([0-9]+)$', label)
        if mat is not None:
            return {
                'group': 'glycan',
                'fragment_type': mat.group(1), 
                'fragment_glycan': mat.group(1) + '0', 
                'loss': 'noloss',
                'charge': int(mat.group(2))
            }
        
        mat = re.match('^([Y])(-[A-Z0-9a-z()]+|\\$)\\+([0-9]+)$', label)
        if mat is not None:
            return {
                'group': 'glycan',
                'fragment_type': mat.group(1), 
                'fragment_glycan': mat.group(1) + mat.group(2), 
                'loss': 'noloss',
                'charge': int(mat.group(3))
            }
        
        return { 
            'group': 'unknown', 
            'fragment_name': label
        }
    
    def parse_matched_ion(self, matched_ion):
        def _parse_matched_ion(s):
            t = s.split('=')
            x = self.parse_ion_annotation(t[0])
            x.update({
                'annotation': t[0]
            })
            if self.matched_ion_has_mz:
                t = t[1].split(',')
                x.update({                    
                    'fragment_mz': float(t[0]),
                    'intensity': float(t[1])
                })
            else:
                x.update({                    
                    'intensity': float(t[1])
                })
            return x
        
        matched_ions_list = [
            _parse_matched_ion(s)
            for s in matched_ion.split(';')
            if len(s) > 0
        ] \
        if matched_ion is not None and \
            isinstance(matched_ion, str) and \
            len(matched_ion) > 0 \
        else []
        
        result = {
            'fragmentAnnotation': [
                x.get('annotation', None) 
                for x in matched_ions_list
            ],
            'fragmentIntensity': [
                x['intensity'] 
                for x in matched_ions_list
            ],
            'fragmentType': [
                x.get('fragment_type', None) 
                for x in matched_ions_list
            ],
            'fragmentNumber': [
                x.get('fragment_number', None) 
                for x in matched_ions_list
            ],
            'fragmentCharge': [
                x.get('charge', None) 
                for x in matched_ions_list
            ],
            'fragmentLossType': [
                x.get('loss', None) 
                for x in matched_ions_list
            ],
            'fragmentGlycan': [
                x.get('fragment_glycan', None) 
                for x in matched_ions_list
            ],
        }
        if self.matched_ion_has_mz:
            result.update({
                'fragmentMZ': [
                    x['fragment_mz'] 
                    for x in matched_ions_list
                ]
            })
        return result
        

    def parse_protein_site(self, protein, protein_site):
        if protein is not None:            
            protein = str(protein).split('/')            
        if protein_site is not None:
            protein_site = str(protein_site).split('/')
            
        if protein is not None and \
            len(protein) > 1 and \
            protein[0] == str(len(protein) - 1):                
            if protein_site is not None and \
                len(protein_site) > 1 and \
                protein_site[0] == protein[0]:
                protein_site.pop(0)
            protein.pop(0)
            
        if protein is not None:
            protein = '/'.join(protein)
        if protein_site is not None:
            protein_site = '/'.join(protein_site)
            
        return protein, protein_site
    
    
    def parse_psm_info(self, psm):
        peptide = psm['Peptide']
        charge = int(psm['Charge'])
        modification = self.parse_modification(psm['Mod'])
        glycan_struct = psm['PlausibleStruct']
        glycan_site = int(psm['GlySite'])
        result = {
            'peptideSequence': peptide,
            'precursorCharge': charge,
            'modification': modification,
            'glycanStruct': glycan_struct,
            'glycanSite': glycan_site
        }

        rt = psm.get('RT', None)
        if rt is not None:
            result.update({
                'rt': float(rt),
            })
        
        precursor_mz = psm.get('PrecursorMZ', None)
        if precursor_mz is not None:
            result.update({
                'precursorMZ': float(precursor_mz),
            })
        
        glycan_fdr = psm.get('GlycanFDR', None)
        peptide_fdr = psm.get('PeptideFDR', None)
        total_fdr = psm.get('TotalFDR', None)
        if total_fdr is not None or \
            peptide_fdr is not None or \
            glycan_fdr is not None:        
            result['metadata'] = {
                'glycanFDR': float(glycan_fdr),
                'peptideFDR': float(peptide_fdr),
                'totalFDR': float(total_fdr)
            }
            
        glycan_score = psm.get('GlyScore', None)
        peptide_score = psm.get('PepScore', None)
        total_score = psm.get('TotalScore', None)
        if total_score is not None or \
            peptide_score is not None or \
            glycan_score is not None:        
            result['metadata'].update({
                'glycanScore': float(glycan_score),
                'peptideScore': float(peptide_score),
                'totalScore': float(total_score)
            })
        
        protein = psm.get('Proteins', None)
        protein_site = psm.get('ProSite', None)
        protein, protein_site = self.parse_protein_site(protein, protein_site)
        if protein is not None:
            result['metadata'].update({
                'protein': protein,
                'proteinSite': protein_site
            })
        
        file = psm.get('RawName', None)
        scan = psm.get('Scan', None)
        if file is not None and scan is not None:
            result['metadata'].update({
                'file': file,
                'scan': int(scan)
            })
        
        return result
    
        
    def psm_to_assay(self, psm):
        result = self.parse_psm_info(psm)        
        fragments = self.parse_matched_ion(psm['MatchedIonInten'])
        result.update({
            'fragments': fragments
        })    
        return result
       
        
    def report_to_assays(self, psm_report, glabel_report, 
                         return_generator=False):
        import pandas as pd
        
        fragment_report = pd.merge(
            psm_report,
            glabel_report \
                .drop(['Peptide'], axis=1) \
                .set_index('Spec'),
            how='inner', 
            left_on=['PepSpec'], 
            right_index=True,
        )

        result = (
            self.psm_to_assay(row) 
            for i, row in fragment_report.iterrows()
        )
        if not return_generator:
            result = list(result)
        return result
        


if __name__ == '__main__':
    converter = pGlycoToAssayConverter()
    
    import pandas as pd
    from io import StringIO
    
    pglyco_str = '''GlySpec	PepSpec	RawName	Scan	RT	PrecursorMH	PrecursorMZ	Charge	Rank	Peptide	Mod	PeptideMH	Glycan(H,N,A,G,F)	PlausibleStruct	GlyID	GlyFrag	GlyMass	GlySite	TotalScore	PepScore	GlyScore	CoreMatched	CoreFuc	MassDeviation	PPM	GlyIonRatio	PepIonRatio	GlyDecoy	PepDecoy	GlycanFDR	PeptideFDR	TotalFDR	Proteins	ProSite
MouseBrain-Z-T-1.9562.9562.2.0.dta	MouseBrain-Z-T-1.9562.9562.2.0.dta	MouseBrain-Z-T-1	9562	2847.60312	3532.40346	1766.70537	2	1	LSVECAJK	5,Carbamidomethyl[C];	920.45058	5 6 1 0 2 	(N(F)(N(H(H(N)(N(H(A))))(H(N(H))(N(F))))))	2058	0 1 0 0 0 ;0 1 0 0 1 ;0 2 0 0 0 ;0 2 0 0 1 ;1 2 0 0 0 ;1 2 0 0 1 ;2 2 0 0 0 ;2 2 0 0 1 ;3 2 0 0 0 ;2 3 0 0 0 ;3 2 0 0 1 ;2 3 0 0 1 ;3 3 0 0 0 ;2 4 0 0 0 ;2 3 0 0 2 ;3 3 0 0 1 ;4 3 0 0 0 ;2 4 0 0 1 ;3 4 0 0 0 ;3 3 1 0 0 ;3 3 0 0 2 ;4 3 0 0 1 ;2 4 0 0 2 ;3 4 0 0 1 ;4 4 0 0 0 ;3 5 0 0 0 ;3 3 1 0 1 ;4 3 1 0 0 ;3 4 1 0 0 ;3 4 0 0 2 ;4 4 0 0 1 ;5 4 0 0 0 ;3 5 0 0 1 ;4 5 0 0 0 ;4 3 1 0 1 ;3 6 0 0 0 ;3 4 1 0 1 ;4 4 1 0 0 ;4 4 0 0 2 ;5 4 0 0 1 ;3 5 0 0 2 ;4 5 0 0 1 ;5 5 0 0 0 ;3 6 0 0 1 ;4 6 0 0 0 ;4 4 1 0 1 ;5 4 1 0 0 ;4 5 1 0 0 ;4 5 0 0 2 ;5 5 0 0 1 ;3 6 0 0 2 ;4 6 0 0 1 ;5 6 0 0 0 ;4 4 1 0 2 ;5 4 1 0 1 ;4 5 1 0 1 ;5 5 1 0 0 ;5 5 0 0 2 ;4 6 1 0 0 ;4 6 0 0 2 ;5 6 0 0 1 ;4 5 1 0 2 ;5 5 1 0 1 ;4 6 1 0 1 ;5 6 1 0 0 ;5 6 0 0 2 ;5 5 1 0 2 ;4 6 1 0 2 ;5 6 1 0 1 ;	2611.95159	7	2.40130	2.17194	2.82724	6	12	0.00129	0.36528	0.05797	0.21429	0	0	0.0478819012157	0.270645491803	0.270645491803	sp|Q62277|SYPH_MOUSE	59
'''

    glabel_str = '''Spec	Peptide	Glycan(H,N,A,G,F)	TheoIon	MatchedIonInten
MouseBrain-Z-T-1.9562.9562.2.0.dta	LSVECAJK	5,6,1,0,2	y7-N1,y7,y6-N1,y6,y5-N1,y5,y4-N1,y4,y3-N1,y3,y2-N1,y2,b7-N1,b7,Y-H5N6A1,Y-H4N4,Y-H3N3F1,Y-N2,Y-H4N5A1,Y-H4N6A1,Y-H3N4F1,Y-H2N3F2,Y-H4N6,Y-H5N6F2,Y-H3N3A1F1,Y-H5N6F1,Y-H2N4F2,Y-H4N3A1F1,Y-H5N5F1,Y-H4N6A1F2,Y-H4N5F1,Y-H5N6,Y-H3N2,Y-H2N4,Y-H4N3,Y-H5N4A1,Y-H5N4F1,Y-H3N5F1,Y-H2N2,Y-H4N5,Y-H3N6F1,Y-H5N5A1F1,Y-H4N4A1F1,Y-H5N5A1,Y-H4N4F1,Y-N1,Y-H5N4,Y-H5N6A1F1,Y-H5N5F2,Y-H2N4F1,Y-H1N2F1,Y-H3N4,Y-H2N3F1,Y-H4N6F1,Y-H3N4A1,Y-H3N3A1,Y-H4N3A1,Y-H3N3,Y-H4N5F2,Y-H5N5,Y-H3N6F2,Y-H4N5A1F1,Y-H3N4F2,Y-H4N5A1F2,Y-H4N4A1F2,Y-H4N4F2,Y-H5N4A1F1,Y-H5N5A1F2,Y-H3N5F2,Y-H2N2F1,Y-H3N3F2,Y-H3N6,Y,Y-H4N4A1,Y-H3N2F1,Y-N2F1,Y-H4N3F1,Y-H3N5,Y-H4N6F2,Y-H1N2,Y-H2N3,Y-H4N6A1F1,Y-H3N4A1F1,Y-N1F1,b1,b2,b3,b4,b5,b6,y1,y$7,y$6,y$5,y$4,y$3,y$2,b$7,Y$,M,	y6+1=857.7;y4+1=2012.5;y2+1=1050.8;Y-N2+1=3452.6;Y-N1+1=21481.3;Y+1=17019.7;Y-N2F1+1=3599.8;Y-N1F1+1=10958.9;Y$+1=5470.2;
'''

    pglyco_report = pd.read_table(StringIO(pglyco_str))
    glabel_report = pd.read_table(StringIO(glabel_str))
    
    print(converter.report_to_assays(
        psm_report=pglyco_report,
        glabel_report=glabel_report
    ))

            